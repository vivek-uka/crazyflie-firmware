/*
 *    ||          ____  _ __
 * +------+      / __ )(_) /_______________ _____  ___
 * | 0xBC |     / __  / / __/ ___/ ___/ __ `/_  / / _ \
 * +------+    / /_/ / / /_/ /__/ /  / /_/ / / /_/  __/
 *  ||  ||    /_____/_/\__/\___/_/   \__,_/ /___/\___/
 *
 * Crazyflie control firmware
 *
 * Copyright (C) 2011-2021 Bitcraze AB
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, in version 3.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 *
 * aideck.c - Deck driver for the AIdeck
 */
#define DEBUG_MODULE "AIDECK"

#include <stdint.h>
#include <string.h>
#include <stdio.h>

#include "stm32fxxx.h"
#include "config.h"
#include "console.h"
#include "uart1.h"
#include "debug.h"
#include "deck.h"
#include "FreeRTOS.h"
#include "task.h"
#include "queue.h"
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "log.h"
#include "param.h"
#include "system.h"
#include "uart1.h"
#include "uart2.h"

static bool isInit = false;
static uint8_t byte;

//Uncomment when NINA printout read is desired from console
//#define DEBUG_NINA_PRINT

#ifdef DEBUG_NINA_PRINT
static void NinaTask(void *param)
{
    systemWaitStart();
    vTaskDelay(M2T(1000));
    DEBUG_PRINT("Starting reading out NINA debugging messages:\n");
    vTaskDelay(M2T(2000));

    // Pull the reset button to get a clean read out of the data
    pinMode(DECK_GPIO_IO4, OUTPUT);
    digitalWrite(DECK_GPIO_IO4, LOW);
    vTaskDelay(10);
    digitalWrite(DECK_GPIO_IO4, HIGH);
    pinMode(DECK_GPIO_IO4, INPUT_PULLUP);

    // Read out the byte the NINA sends and immediately send it to the console.
    uint8_t byte;
    while (1)
    {
        if (uart2GetDataWithDefaultTimeout(&byte) == true)
        {
            consolePutchar(byte);
        }
    }
}
#endif

static void Gap8Task(void *param)
{
    systemWaitStart();
    vTaskDelay(M2T(1000));

    // Pull the reset button to get a clean read out of the data
    // pinMode(DECK_GPIO_IO4, OUTPUT);
    // digitalWrite(DECK_GPIO_IO4, LOW);
    // vTaskDelay(10);
    // digitalWrite(DECK_GPIO_IO4, HIGH);
    // pinMode(DECK_GPIO_IO4, INPUT_PULLUP);

    // static logVarId_t logIdacc_x;
    // static logVarId_t logIdacc_y;
    static logVarId_t logHeight;
    static logVarId_t logIdacc_z;
    static logVarId_t logIdgyro_x;
    static logVarId_t logIdgyro_y;
    static logVarId_t logIdgyro_z;

    // get the values from log framework
    // logIdacc_x = logGetVarId("acc", "x");
    // logIdacc_y = logGetVarId("acc", "y");

    logHeight = logGetVarId("kalman", "stateZ");
    logIdacc_z = logGetVarId("acc", "z");
    logIdgyro_x = logGetVarId("gyro", "x");
    logIdgyro_y = logGetVarId("gyro", "y");
    logIdgyro_z = logGetVarId("gyro", "z");

    unsigned char packet[28];
    packet[0] = 0xfe;
    packet[1] = 0xfe;
    packet[26] = 0xfe;
    packet[27] = 0xfe;

    while (1)
    {
        // Read out the byte the Gap8 sends and immediately send it to the console.
        // uart1GetDataWithDefaultTimeout(&byte);

        // float acc_x = logGetFloat(logIdacc_x);
        // float acc_y = logGetFloat(logIdacc_y);

        float height =  logGetFloat(logHeight);

        float acc_z = logGetFloat(logIdacc_z);
        float gyro_x = logGetFloat(logIdgyro_x);
        float gyro_y = logGetFloat(logIdgyro_y);
        float gyro_z = logGetFloat(logIdgyro_z);

        // get current os time tick
        uint32_t osTick = xTaskGetTickCount();   // 32 bits --> 4 bytes -- the same with "float"
        float time_tick = (float) osTick;

        union {
            float a;
            unsigned char bytes[4];
        } thing;
        thing.a = time_tick;                     // test

        for(uint16_t i=0; i<4; i++)
        {
            packet[2+i] = thing.bytes[i];
        }
        thing.a = height;
        for(uint16_t i=0; i<4; i++)
        {
            packet[6+i] = thing.bytes[i];
        }
        thing.a = acc_z;
        for(uint16_t i=0; i<4; i++)
        {
            packet[10+i] = thing.bytes[i];
        }
        thing.a = gyro_x;
        for(uint16_t i=0; i<4; i++)
        {
            packet[14+i] = thing.bytes[i];
        }
        thing.a = gyro_y;
        for(uint16_t i=0; i<4; i++)
        {
            packet[18+i] = thing.bytes[i];
        }
        thing.a = gyro_z;
        for(uint16_t i=0; i<4; i++)
        {
            packet[22+i] = thing.bytes[i];
        }
        uart1SendData(28, packet);
        // vTaskDelay(M2T((rand()%10) + 20));
        vTaskDelay(M2T(40));
    }
}

static void aideckInit(DeckInfo *info)
{

    if (isInit)
        return;

    // Intialize the UART for the GAP8
    uart1Init(115200);
    // Initialize task for the GAP8
    xTaskCreate(Gap8Task, AI_DECK_GAP_TASK_NAME, AI_DECK_TASK_STACKSIZE, NULL,
                AI_DECK_TASK_PRI, NULL);

#ifdef DEBUG_NINA_PRINT
    // Initialize the UART for the NINA
    uart2Init(115200);
    // Initialize task for the NINA
    xTaskCreate(NinaTask, AI_DECK_NINA_TASK_NAME, AI_DECK_TASK_STACKSIZE, NULL,
                AI_DECK_TASK_PRI, NULL);

#endif

    isInit = true;
}

static bool aideckTest()
{

    return true;
}

static const DeckDriver aideck_deck = {
    .vid = 0xBC,
    .pid = 0x12,
    .name = "bcAI",

    .usedGpio = DECK_USING_IO_4,
    .usedPeriph = DECK_USING_UART1,

    .init = aideckInit,
    .test = aideckTest,
};

LOG_GROUP_START(aideck)
LOG_ADD(LOG_UINT8, receivebyte, &byte)
LOG_GROUP_STOP(aideck)

/** @addtogroup deck
*/
PARAM_GROUP_START(deck)

/**
 * @brief Nonzero if [AI deck](%https://store.bitcraze.io/collections/decks/products/ai-deck-1-1) is attached
 */
PARAM_ADD_CORE(PARAM_UINT8 | PARAM_RONLY, bcAIDeck, &isInit)

PARAM_GROUP_STOP(deck)

DECK_DRIVER(aideck_deck);
