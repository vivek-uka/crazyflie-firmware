/*
 *    ||          ____  _ __
 * +------+      / __ )(_) /_______________ _____  ___
 * | 0xBC |     / __  / / __/ ___/ ___/ __ `/_  / / _ \
 * +------+    / /_/ / / /_/ /__/ /  / /_/ / / /_/  __/
 *  ||  ||    /_____/_/\__/\___/_/   \__,_/ /___/\___/
 *
 * Crazyflie firmware.
 *
 * Copyright 2018, Bitcraze AB
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * tdoaEngine.c is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with tdoaEngine.c. If not, see <http://www.gnu.org/licenses/>.
 */


/*
Implementation of LPS TDoA Tag functionality

The tag is assumed to move around in a large system of anchors. Any anchor ids
can be used, and the same anchor id can even be used by multiple anchors as long
as they are not visible in the same area. It is assumed that the anchor density
is evenly distributed in the covered volume and that 5-20 anchors are visible
in every point. The tag is attached to a physical object and the expected
velocity is a few m/s, this means that anchors are within range for a time
period of seconds.

The implementation must handle
1. An infinite number of anchors, where around 20 are visible at one time
2. Any anchor ids
3. Dynamically changing visibility of anchors over time
4. Random TX times from anchors with possible packet collisions and packet loss

*/

#include <string.h>

#define DEBUG_MODULE "TDOA_ENGINE"
#include "debug.h"

#include "tdoaEngine.h"
#include "tdoaStats.h"
#include "clockCorrectionEngine.h"
#include "physicalConstants.h"

#include "log.h"

#define MEASUREMENT_NOISE_STD 0.15f

// define the locodeck_freq for easy of use in this code
// it is the same with the "#define LOCODECK_TS_FREQ (499.2e6 * 128)" in locodeck.h
#define LOCODECK_FREQ (499.2e6 * 128)
#define ANTENNA_OFFSET 154.6
#define LOCODECK_ANTENNA_DELAY  ((ANTENNA_OFFSET * LOCODECK_FREQ) / SPEED_OF_LIGHT) // In radio ticks
// log data
static float log_tdoa_in_clT = 0.0f;
static float log_d1 = 0.0f;
static float log_d2 = 0.0f;
static uint16_t log_tof = 0;
static float tof_dist = 0.0f;

void tdoaEngineInit(tdoaEngineState_t* engineState, const uint32_t now_ms, tdoaEngineSendTdoaToEstimator sendTdoaToEstimator, const double locodeckTsFreq, const tdoaEngineMatchingAlgorithm_t matchingAlgorithm) {
  tdoaStorageInitialize(engineState->anchorInfoArray);
  tdoaStatsInit(&engineState->stats, now_ms);
  engineState->sendTdoaToEstimator = sendTdoaToEstimator;   // for TDOA3, init the function to sendTdoaToEstimatorCallback in lpsTdoa3Tag.c
  engineState->locodeckTsFreq = locodeckTsFreq;
  engineState->matchingAlgorithm = matchingAlgorithm;

  engineState->matching.offset = 0;
}

static void enqueueTDOA(const tdoaAnchorContext_t* anchorACtx, const tdoaAnchorContext_t* anchorBCtx, double distanceDiff, tdoaEngineState_t* engineState) {
  tdoaStats_t* stats = &engineState->stats;
  
  tdoaMeasurement_t tdoa = {
    .stdDev = MEASUREMENT_NOISE_STD,
    .distanceDiff = distanceDiff
  };

  if (tdoaStorageGetAnchorPosition(anchorACtx, &tdoa.anchorPositions[0]) && tdoaStorageGetAnchorPosition(anchorBCtx, &tdoa.anchorPositions[1])) {
    STATS_CNT_RATE_EVENT(&stats->packetsToEstimator);

    uint8_t idA = tdoaStorageGetId(anchorACtx);
    uint8_t idB = tdoaStorageGetId(anchorBCtx);
    if (idA == stats->anchorId && idB == stats->remoteAnchorId) {
      stats->tdoa = distanceDiff;
    }
    if (idB == stats->anchorId && idA == stats->remoteAnchorId) {
      stats->tdoa = -distanceDiff;
    }
    tdoa.anchorIds[0] = idA;
    tdoa.anchorIds[1] = idB;
    engineState->sendTdoaToEstimator(&tdoa);   // for TDOA3, call the function sendTdoaToEstimatorCallback
  }
}

static bool updateClockCorrection(tdoaAnchorContext_t* anchorCtx, const int64_t txAn_in_cl_An, const int64_t rxAn_by_T_in_cl_T, tdoaStats_t* stats) {
  bool sampleIsReliable = false;

  const int64_t latest_rxAn_by_T_in_cl_T = tdoaStorageGetRxTime(anchorCtx);
  const int64_t latest_txAn_in_cl_An = tdoaStorageGetTxTime(anchorCtx);

  if (latest_rxAn_by_T_in_cl_T != 0 && latest_txAn_in_cl_An != 0) {
    double clockCorrectionCandidate = clockCorrectionEngineCalculate(rxAn_by_T_in_cl_T, latest_rxAn_by_T_in_cl_T, txAn_in_cl_An, latest_txAn_in_cl_An, TDOA_ENGINE_TRUNCATE_TO_ANCHOR_TS_BITMAP);
    sampleIsReliable = clockCorrectionEngineUpdate(tdoaStorageGetClockCorrectionStorage(anchorCtx), clockCorrectionCandidate);

    if (sampleIsReliable){
      if (tdoaStorageGetId(anchorCtx) == stats->anchorId) {
        stats->clockCorrection = tdoaStorageGetClockCorrection(anchorCtx);
        STATS_CNT_RATE_EVENT(&stats->clockCorrectionCount);
      }
    }
  }

  return sampleIsReliable;
}

// [change] log every value to compute the TDOA meas.
// We do not know the start time of tag and anchor, so that we can not compute two TOA separately and then get the tdoa meas.
static int64_t calcTDoA(const tdoaAnchorContext_t* otherAnchorCtx, const tdoaAnchorContext_t* anchorCtx, const int64_t txAn_in_cl_An, const int64_t rxAn_by_T_in_cl_T) {
  const uint8_t otherAnchorId = tdoaStorageGetId(otherAnchorCtx);

  const int64_t tof_Ar_to_An_in_cl_An = tdoaStorageGetTimeOfFlight(anchorCtx, otherAnchorId);
  const int64_t rxAr_by_An_in_cl_An = tdoaStorageGetRemoteRxTime(anchorCtx, otherAnchorId);
  const double clockCorrection = tdoaStorageGetClockCorrection(anchorCtx);

  const int64_t rxAr_by_T_in_cl_T = tdoaStorageGetRxTime(otherAnchorCtx);

  const int64_t delta_txAr_to_txAn_in_cl_An = (tof_Ar_to_An_in_cl_An + tdoaEngineTruncateToAnchorTimeStamp(txAn_in_cl_An - rxAr_by_An_in_cl_An));
  const int64_t timeDiffOfArrival_in_cl_T =  tdoaEngineTruncateToAnchorTimeStamp(rxAn_by_T_in_cl_T - rxAr_by_T_in_cl_T) - delta_txAr_to_txAn_in_cl_An  * clockCorrection;

  // log data
  // otherAnchorCtx is anA,  anchorCtx is anB
  uint8_t idA = tdoaStorageGetId(otherAnchorCtx);
  uint8_t idB = tdoaStorageGetId(anchorCtx);
  if(idA==(uint8_t)1 && idB == (uint8_t)2){
      // d1_in_cl_T = t1 - alpha * t_an_r + alpha * tof
    //   int64_t d1_in_cl_T = tdoaEngineTruncateToAnchorTimeStamp(rxAr_by_T_in_cl_T) - clockCorrection*tdoaEngineTruncateToAnchorTimeStamp(rxAr_by_An_in_cl_An) + clockCorrection * tof_Ar_to_An_in_cl_An;
      // d2_in_cl_T = t2 - alpha * t_an_t  
    //   int64_t d2_in_cl_T = tdoaEngineTruncateToAnchorTimeStamp(rxAn_by_T_in_cl_T) - clockCorrection*tdoaEngineTruncateToAnchorTimeStamp(txAn_in_cl_An);

      // the values in dummy_d1 and dummy_d2 are correct    
    //   int64_t dummy_d1 = tdoaEngineTruncateToAnchorTimeStamp(rxAn_by_T_in_cl_T - rxAr_by_T_in_cl_T);
    //   int64_t dummy_d2 = delta_txAr_to_txAn_in_cl_An  * clockCorrection;
      log_tof = (uint16_t)tof_Ar_to_An_in_cl_An;
      // need to minus the ANTENNA_OFFSET (see tdoa3_tof.py in anchor sniffer code)
      // now tof_dist is correct
      tof_dist = SPEED_OF_LIGHT * log_tof / LOCODECK_FREQ - ANTENNA_OFFSET;

      // the values in dmy_d1 and dmy_d2 are correct 
      // when I use the dmy_d1 to compute distance, the value has a large base (base + distance)
      // yte compute the distance_diff is correct
      int64_t dmy_d1 = tdoaEngineTruncateToAnchorTimeStamp(rxAr_by_T_in_cl_T) - clockCorrection*tdoaEngineTruncateToAnchorTimeStamp(rxAr_by_An_in_cl_An) + clockCorrection * tof_Ar_to_An_in_cl_An;
      int64_t dmy_d2 = tdoaEngineTruncateToAnchorTimeStamp(rxAn_by_T_in_cl_T) - clockCorrection*tdoaEngineTruncateToAnchorTimeStamp(txAn_in_cl_An);
      
      // when put into the logging (float), we lost some values. (for both dmy_d1 and dummy_d1)
      log_d1 = dmy_d1;
      log_d2 = dmy_d2;
    //   log_d1 = dummy_d1;
    //   log_d2 = dummy_d2;
      log_tdoa_in_clT  = SPEED_OF_LIGHT * timeDiffOfArrival_in_cl_T / LOCODECK_FREQ;
  } 
  return timeDiffOfArrival_in_cl_T;
}

// [change] log every value to compute the TDOA meas.
// locodeckTsFreq is init in locodeck.h. ----->  #define LOCODECK_TS_FREQ (499.2e6 * 128)
static double calcDistanceDiff(const tdoaAnchorContext_t* otherAnchorCtx, const tdoaAnchorContext_t* anchorCtx, const int64_t txAn_in_cl_An, const int64_t rxAn_by_T_in_cl_T, const double locodeckTsFreq) {
  const int64_t tdoa = calcTDoA(otherAnchorCtx, anchorCtx, txAn_in_cl_An, rxAn_by_T_in_cl_T);
  return SPEED_OF_LIGHT * tdoa / locodeckTsFreq;
}

static bool matchRandomAnchor(tdoaEngineState_t* engineState, tdoaAnchorContext_t* otherAnchorCtx, const tdoaAnchorContext_t* anchorCtx, const bool doExcludeId, const uint8_t excludedId) {
  engineState->matching.offset++;
  int remoteCount = 0;
  tdoaStorageGetRemoteSeqNrList(anchorCtx, &remoteCount, engineState->matching.seqNr, engineState->matching.id);

  uint32_t now_ms = anchorCtx->currentTime_ms;

  // Loop over the candidates and pick the first one that is useful
  // An offset (updated for each call) is added to make sure we start at
  // different positions in the list and vary which candidate to choose
  for (int i = engineState->matching.offset; i < (remoteCount + engineState->matching.offset); i++) {
    uint8_t index = i % remoteCount;
    const uint8_t candidateAnchorId = engineState->matching.id[index];
    if (!doExcludeId || (excludedId != candidateAnchorId)) {
      if (tdoaStorageGetCreateAnchorCtx(engineState->anchorInfoArray, candidateAnchorId, now_ms, otherAnchorCtx)) {
        if (engineState->matching.seqNr[index] == tdoaStorageGetSeqNr(otherAnchorCtx) && tdoaStorageGetTimeOfFlight(anchorCtx, candidateAnchorId)) {
          return true;
        }
      }
    }
  }

  otherAnchorCtx->anchorInfo = 0;
  return false;
}

static bool matchYoungestAnchor(tdoaEngineState_t* engineState, tdoaAnchorContext_t* otherAnchorCtx, const tdoaAnchorContext_t* anchorCtx, const bool doExcludeId, const uint8_t excludedId) {
    int remoteCount = 0;
    tdoaStorageGetRemoteSeqNrList(anchorCtx, &remoteCount, engineState->matching.seqNr, engineState->matching.id);

    uint32_t now_ms = anchorCtx->currentTime_ms;
    uint32_t youmgestUpdateTime = 0;
    int bestId = -1;

    for (int index = 0; index < remoteCount; index++) {
      const uint8_t candidateAnchorId = engineState->matching.id[index];
      if (!doExcludeId || (excludedId != candidateAnchorId)) {
        if (tdoaStorageGetTimeOfFlight(anchorCtx, candidateAnchorId)) {
          if (tdoaStorageGetCreateAnchorCtx(engineState->anchorInfoArray, candidateAnchorId, now_ms, otherAnchorCtx)) {
            uint32_t updateTime = otherAnchorCtx->anchorInfo->lastUpdateTime;
            if (updateTime > youmgestUpdateTime) {
              if (engineState->matching.seqNr[index] == tdoaStorageGetSeqNr(otherAnchorCtx)) {
                youmgestUpdateTime = updateTime;
                bestId = candidateAnchorId;
              }
            }
          }
        }
      }
    }

    if (bestId >= 0) {
      tdoaStorageGetCreateAnchorCtx(engineState->anchorInfoArray, bestId, now_ms, otherAnchorCtx);
      return true;
    }

    otherAnchorCtx->anchorInfo = 0;
    return false;
}

static bool findSuitableAnchor(tdoaEngineState_t* engineState, tdoaAnchorContext_t* otherAnchorCtx, const tdoaAnchorContext_t* anchorCtx, const bool doExcludeId, const uint8_t excludedId) {
  bool result = false;

  if (tdoaStorageGetClockCorrection(anchorCtx) > 0.0) {
    switch(engineState->matchingAlgorithm) {
      case TdoaEngineMatchingAlgorithmRandom:
        result = matchRandomAnchor(engineState, otherAnchorCtx, anchorCtx, doExcludeId, excludedId);
        break;

      case TdoaEngineMatchingAlgorithmYoungest:
        result = matchYoungestAnchor(engineState, otherAnchorCtx, anchorCtx, doExcludeId, excludedId);
        break;

      default:
        // Do nothing
        break;
    }
  }

  return result;
}

void tdoaEngineGetAnchorCtxForPacketProcessing(tdoaEngineState_t* engineState, const uint8_t anchorId, const uint32_t currentTime_ms, tdoaAnchorContext_t* anchorCtx) {
  if (tdoaStorageGetCreateAnchorCtx(engineState->anchorInfoArray, anchorId, currentTime_ms, anchorCtx)) {
    STATS_CNT_RATE_EVENT(&engineState->stats.contextHitCount);
  } else {
    STATS_CNT_RATE_EVENT(&engineState->stats.contextMissCount);
  }
}

void tdoaEngineProcessPacket(tdoaEngineState_t* engineState, tdoaAnchorContext_t* anchorCtx, const int64_t txAn_in_cl_An, const int64_t rxAn_by_T_in_cl_T) {
  tdoaEngineProcessPacketFiltered(engineState, anchorCtx, txAn_in_cl_An, rxAn_by_T_in_cl_T, false, 0);
}

void tdoaEngineProcessPacketFiltered(tdoaEngineState_t* engineState, tdoaAnchorContext_t* anchorCtx, const int64_t txAn_in_cl_An, const int64_t rxAn_by_T_in_cl_T, const bool doExcludeId, const uint8_t excludedId) {
  bool timeIsGood = updateClockCorrection(anchorCtx, txAn_in_cl_An, rxAn_by_T_in_cl_T, &engineState->stats);
  if (timeIsGood) {
    STATS_CNT_RATE_EVENT(&engineState->stats.timeIsGood);

    tdoaAnchorContext_t otherAnchorCtx;
    if (findSuitableAnchor(engineState, &otherAnchorCtx, anchorCtx, doExcludeId, excludedId)) {
      STATS_CNT_RATE_EVENT(&engineState->stats.suitableDataFound);
      double tdoaDistDiff = calcDistanceDiff(&otherAnchorCtx, anchorCtx, txAn_in_cl_An, rxAn_by_T_in_cl_T, engineState->locodeckTsFreq);
      // otherAnchorCtx is anA,  anchorCtx is anB
      enqueueTDOA(&otherAnchorCtx, anchorCtx, tdoaDistDiff, engineState);
    }
  }
}

LOG_GROUP_START(engine)
LOG_ADD(LOG_FLOAT,  tdoa_in_clT, &log_tdoa_in_clT)
LOG_ADD(LOG_FLOAT,  d1_in_clT,   &log_d1)
LOG_ADD(LOG_FLOAT,  d2_in_clT,   &log_d2)
LOG_ADD(LOG_UINT16, tof_in_cl,   &log_tof)
LOG_ADD(LOG_FLOAT,  dist-tof,    &tof_dist)
LOG_GROUP_STOP(engine)

