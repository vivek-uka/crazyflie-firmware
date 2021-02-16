// DSL title

#pragma once

#include "kalman_core.h"

// M-estimation based robust Kalman filter update for UWB TWR measurements
void kalmanCoreRobustUpdateWithDistance(kalmanCoreData_t* this, distanceMeasurement_t *d);
