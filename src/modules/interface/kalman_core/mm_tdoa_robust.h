// DSL title

#pragma once

#include "kalman_core.h"

// M-estimation based robust Kalman filter update for UWB TDOA measurements
void kalmanCoreRobustUpdateWithTDOA(kalmanCoreData_t* this, tdoaMeasurement_t *tdoa);



