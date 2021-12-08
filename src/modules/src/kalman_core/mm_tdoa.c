/**
 * ,---------,       ____  _ __
 * |  ,-^-,  |      / __ )(_) /_______________ _____  ___
 * | (  O  ) |     / __  / / __/ ___/ ___/ __ `/_  / / _ \
 * | / ,--'  |    / /_/ / / /_/ /__/ /  / /_/ / / /_/  __/
 *    +------`   /_____/_/\__/\___/_/   \__,_/ /___/\___/
 *
 * Crazyflie control firmware
 *
 * Copyright (C) 2021 Bitcraze AB
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
 */

#include "mm_tdoa.h"
#include "outlierFilter.h"
#include "test_support.h"

// TODO krri What is this used for? Do we still need it?
TESTABLE_STATIC uint32_t tdoaCount = 0;

void kalmanCoreUpdateWithTDOA(kalmanCoreData_t* this, tdoaMeasurement_t *tdoa)
{
  if (tdoaCount >= 100)
  {
    /**
     * Measurement equation:
     * dR = dT + d1 - d0
     */

    float measurement = tdoa->distanceDiff;
    // predict based on current state
    float x = this->S[KC_STATE_X];
    float y = this->S[KC_STATE_Y];
    float z = this->S[KC_STATE_Z];

    // consider lever-arm effect for CF_Bolt
    float t_uv[3] = {-0.01245, 0.00127, 0.0908};     // translation vector from the quadrotor to UWB tag
    // compute the position of uwb tag: p_uwb = R_iv.dot(t_uv) + Xpr
    float p_uwb[3] = {0}; 
    p_uwb[0] = this->R[0][0] * t_uv[0] + this->R[0][1] * t_uv[1] + this->R[0][2] * t_uv[2] + x;
    p_uwb[1] = this->R[1][0] * t_uv[0] + this->R[1][1] * t_uv[1] + this->R[1][2] * t_uv[2] + y;
    p_uwb[2] = this->R[2][0] * t_uv[0] + this->R[2][1] * t_uv[1] + this->R[2][2] * t_uv[2] + z;

    // anchor positions
    float x1 = tdoa->anchorPositions[1].x, y1 = tdoa->anchorPositions[1].y, z1 = tdoa->anchorPositions[1].z;
    float x0 = tdoa->anchorPositions[0].x, y0 = tdoa->anchorPositions[0].y, z0 = tdoa->anchorPositions[0].z;
    
    float dx1 = p_uwb[0] - x1;
    float dy1 = p_uwb[1] - y1;
    float dz1 = p_uwb[2] - z1;

    float dx0 = p_uwb[0] - x0;
    float dy0 = p_uwb[1] - y0;
    float dz0 = p_uwb[2] - z0;

    float d1 = sqrtf(powf(dx1, 2) + powf(dy1, 2) + powf(dz1, 2));
    float d0 = sqrtf(powf(dx0, 2) + powf(dy0, 2) + powf(dz0, 2));

    float predicted = d1 - d0;
    float error = measurement - predicted;

    float h[KC_STATE_DIM] = {0};
    arm_matrix_instance_f32 H = {1, KC_STATE_DIM, h};

    if ((d0 != 0.0f) && (d1 != 0.0f)) {
      // g_p
      h[KC_STATE_X] = (dx1 / d1 - dx0 / d0);
      h[KC_STATE_Y] = (dy1 / d1 - dy0 / d0);
      h[KC_STATE_Z] = (dz1 / d1 - dz0 / d0);
      // g_v ={0,0,0}
      // compute g_q
      float d_vec[3] = {0};
      d_vec[0] = this->q[0]*t_uv[0] - this->q[3]*t_uv[1] + this->q[2]*t_uv[2];
      d_vec[1] = this->q[0]*t_uv[1] + this->q[3]*t_uv[0] - this->q[1]*t_uv[2];
      d_vec[2] = this->q[0]*t_uv[2] - this->q[2]*t_uv[0] + this->q[1]*t_uv[1];

      // q_v.T.dot(t_uv)
      float qt = this->q[1]*t_uv[0] + this->q[2]*t_uv[1] + this->q[3]*t_uv[2];
      float d_mat[3][3] = {0};
      d_mat[0][0] = qt + this->q[1]*t_uv[0] - t_uv[0]*this->q[1];
      d_mat[0][1] = this->q[1]*t_uv[1] - t_uv[0]*this->q[2] + this->q[0] * t_uv[2];
      d_mat[0][2] = this->q[1]*t_uv[2] - t_uv[0]*this->q[3] - this->q[0] * t_uv[1];

      d_mat[1][0] = this->q[2]*t_uv[0] - t_uv[1]*this->q[1] - this->q[0] * t_uv[2];
      d_mat[1][1] = qt + this->q[2]*t_uv[1] - t_uv[1]*this->q[2];
      d_mat[1][2] = this->q[2]*t_uv[2] - t_uv[1]*this->q[3] + this->q[0] * t_uv[0];

      d_mat[2][0] = this->q[3]*t_uv[0] - t_uv[2]*this->q[1] + this->q[0] * t_uv[1];
      d_mat[2][1] = this->q[3]*t_uv[1] - t_uv[2]*this->q[2] - this->q[0] * t_uv[0];
      d_mat[2][2] = qt + this->q[3]*t_uv[2] - t_uv[2]*this->q[3];
      // compute d_RVq

      // compute H_x

      // compute H_dx

      // compute H

      kalmanCoreScalarUpdate(this, &H, error, tdoa->stdDev);

      // very coarse outlier detection
      // vector_t jacobian = {
      //   .x = h[KC_STATE_X],
      //   .y = h[KC_STATE_Y],
      //   .z = h[KC_STATE_Z],
      // };

      // point_t estimatedPosition = {
      //   .x = this->S[KC_STATE_X],
      //   .y = this->S[KC_STATE_Y],
      //   .z = this->S[KC_STATE_Z],
      // };

      // bool sampleIsGood = outlierFilterValidateTdoaSteps(tdoa, error, &jacobian, &estimatedPosition);
      // if (sampleIsGood) {
      //   kalmanCoreScalarUpdate(this, &H, error, tdoa->stdDev);
      // }

    }
  }

  tdoaCount++;
}
