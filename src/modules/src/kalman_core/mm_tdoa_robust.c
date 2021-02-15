// DSL title

#include "mm_tdoa_robust.h"
#include "test_support.h"
     

#define MAX_ITER (2) // maximum iteration is set to 2. 

// Cholesky Decomposition for a nxn psd matrix(from scratch)
// Reference: https://www.geeksforgeeks.org/cholesky-decomposition-matrix-decomposition/
static void Cholesky_Decomposition(int n, float matrix[n][n],  float lower[n][n]){
    // Decomposing a matrix into Lower Triangular 
    for (int i = 0; i < n; i++) { 
        for (int j = 0; j <= i; j++) { 
            float sum = 0.0; 
            if (j == i) // summation for diagnols 
            { 
                for (int k = 0; k < j; k++) 
                    sum += powf(lower[j][k], 2); 
                lower[j][j] = sqrtf(matrix[j][j] - sum); 
            } else { 
                // Evaluating L(i, j) using L(j, j) 
                for (int k = 0; k < j; k++) 
                    sum += (lower[i][k] * lower[j][k]); 
                lower[i][j] = (matrix[i][j] - sum) / lower[j][j]; 
            } 
        } 
    }
} 

// Weight function for GM Robust cost function
static void GM_UWB(float e, float * GM_e){
    float sigma = 2.0;                        // param1: 2.0,   param2: 2.5
    float GM_dn = sigma + e*e;
    *GM_e = (sigma * sigma)/(GM_dn * GM_dn);
}

static void GM_state(float e, float * GM_e){
    float sigma = 1.5;                        // param1: 1.5
    float GM_dn = sigma + e*e;
    *GM_e = (sigma * sigma)/(GM_dn * GM_dn);
}

TESTABLE_STATIC uint32_t tdoaCount = 0;
// robsut update function
void kalmanCoreRobustUpdateWithTDOA(kalmanCoreData_t* this, tdoaMeasurement_t *tdoa)
{
    if (tdoaCount >= 100)
  {
    /**
     * Measurement equation:
     * dR = dT + d1 - d0
     */
	float measurement = 0.0f;
    // x prior (x_check)
    float x = this->S[KC_STATE_X];   float y = this->S[KC_STATE_Y];   float z = this->S[KC_STATE_Z];

    float x1 = tdoa->anchorPosition[1].x, y1 = tdoa->anchorPosition[1].y, z1 = tdoa->anchorPosition[1].z;
    float x0 = tdoa->anchorPosition[0].x, y0 = tdoa->anchorPosition[0].y, z0 = tdoa->anchorPosition[0].z;

    float dx1 = x - x1;     float dy1 = y - y1;     float dz1 = z - z1;
    float dx0 = x - x0;     float dy0 = y - y0;     float dz0 = z - z0;

    float d1 = sqrtf(powf(dx1, 2) + powf(dy1, 2) + powf(dz1, 2));
    float d0 = sqrtf(powf(dx0, 2) + powf(dy0, 2) + powf(dz0, 2));
    // if measurements make sense 
    if ((d0 != 0.0f) && (d1 != 0.0f)) {
        
        float predicted = d1 - d0;
        // without bias compensation
        measurement = tdoa->distanceDiff;
        
        // innovation term based on x_check
        float error_check = measurement - predicted;    // error_check

        // ------------------------ unit test ------------------------ //
        // measurement = 2.45f;
        // error_check = 0.12f;
        // S[0]=1.5f;   S[1]=0.0f; S[2]=0.0f; S[3]=0.0f; S[4]=0.0f; S[5]=0.0f; S[6]=0.0f; S[7]=0.0f; S[8]=0.0f;
        // static float P[9][9]={0};
        // P[0][0] = 1.0f;  P[1][1] = 1.0f;  P[2][2] = 1.0f;
        // P[3][3] = 0.1f;  P[4][4] = 0.1f;  P[5][5] = 0.1f;
        // P[6][6] = 0.1f;  P[7][7] = 0.1f;  P[8][8] = 0.1f;
        // ---------------------- matrix defination ----------------------------- //
        static float P_chol[KC_STATE_DIM][KC_STATE_DIM]; 
        static arm_matrix_instance_f32 Pc_m = {KC_STATE_DIM, KC_STATE_DIM, (float *)P_chol};
        // Pc.T
        static float Pc_tran[KC_STATE_DIM][KC_STATE_DIM];
        static arm_matrix_instance_f32 Pc_tran_m = {KC_STATE_DIM, KC_STATE_DIM, (float *)Pc_tran};

        float h[KC_STATE_DIM] = {0};
        arm_matrix_instance_f32 H = {1, KC_STATE_DIM, h};    // H is a row vector

        static float Kw[KC_STATE_DIM];           // The Kalman gain as a column vector
        static arm_matrix_instance_f32 Kwm = {KC_STATE_DIM, 1, (float *)Kw};

        // float error_x[STATE_DIM]={0}; // is not useful
        // arm_matrix_instance_f32 error_x_mat = {STATE_DIM, 1, error_x};
        
        static float e_x[KC_STATE_DIM];
        static arm_matrix_instance_f32 e_x_m = {KC_STATE_DIM, 1, e_x};
        
        static float Pc_inv[KC_STATE_DIM][KC_STATE_DIM];
        static arm_matrix_instance_f32 Pc_inv_m = {KC_STATE_DIM, KC_STATE_DIM, (float *)Pc_inv};
        
        // rescale matrix
        static float wx_inv[KC_STATE_DIM][KC_STATE_DIM];
        static arm_matrix_instance_f32 wx_invm = {KC_STATE_DIM, KC_STATE_DIM, (float *)wx_inv};
        // tmp matrix for P_chol inverse
        static float tmp1[KC_STATE_DIM][KC_STATE_DIM];
        static arm_matrix_instance_f32 tmp1m = {KC_STATE_DIM, KC_STATE_DIM, (float *)tmp1};

        static float Pc_w_inv[KC_STATE_DIM][KC_STATE_DIM];
        static arm_matrix_instance_f32 Pc_w_invm = {KC_STATE_DIM, KC_STATE_DIM, (float *)Pc_w_inv};

        static float P_w[KC_STATE_DIM][KC_STATE_DIM];
        static arm_matrix_instance_f32 P_w_m = {KC_STATE_DIM, KC_STATE_DIM, (float *)P_w};

        static float HTd[KC_STATE_DIM];
        static arm_matrix_instance_f32 HTm = {KC_STATE_DIM, 1, HTd};

        static float PHTd[KC_STATE_DIM];
        static arm_matrix_instance_f32 PHTm = {KC_STATE_DIM, 1, PHTd};
        // ------------------- Some initialization -----------------------//
        // float xpr[STATE_DIM] = {0.0};                // x prior (error state), set to be zeros 
        static float x_err[KC_STATE_DIM] = {0.0};          // x_err comes from the KF update is the state of error state Kalman filter, set to be zero initially
        static arm_matrix_instance_f32 x_errm = {KC_STATE_DIM, 1, x_err};
        static float X_state[KC_STATE_DIM] = {0.0};
        float P_iter[KC_STATE_DIM][KC_STATE_DIM];
        matrixcopy(KC_STATE_DIM, KC_STATE_DIM, P_iter,this->P);     // init P_iter as P_prior
        // construct Q
        float Q_iter = tdoa->stdDev * tdoa->stdDev;
        vectorcopy(KC_STATE_DIM, X_state, this->S);             // copy Xpr to X_State and then update in each iterations

        // ---------------------- Start iteration ----------------------- //
        for (int iter = 0; iter < MAX_ITER; iter++){
            // apply cholesky decomposition for the prior covariance matrix 
            Cholesky_Decomposition(KC_STATE_DIM, P_iter, P_chol);               // P_chol is a lower triangular matrix
            mat_trans(&Pc_m, &Pc_tran_m);
            // test cholesky --> cholesky is correct !!!!
            //----------------------------------------------//
            // decomposition for measurement covariance (scalar case)
            // only true in scalar case
            float Q_chol = sqrtf(Q_iter);           
            // construct H matrix
            // X_state updates in each iteration
            float x_iter = X_state[KC_STATE_X];   float y_iter = X_state[KC_STATE_Y];   float z_iter = X_state[KC_STATE_Z];

            float x1 = tdoa->anchorPosition[1].x, y1 = tdoa->anchorPosition[1].y, z1 = tdoa->anchorPosition[1].z;
            float x0 = tdoa->anchorPosition[0].x, y0 = tdoa->anchorPosition[0].y, z0 = tdoa->anchorPosition[0].z;     

            float dx1 = x_iter - x1;  float dy1 = y_iter - y1;  float dz1 = z_iter - z1;
            float dx0 = x_iter - x0;  float dy0 = y_iter - y0;  float dz0 = z_iter - z0;
            // ---------------- debug inputs ---------------- //
            // dx1 = x_iter + 3.0f;  dy1 = y_iter - 1.0f;  dz1 = z_iter - 2.0f;
            // dx0 = x_iter - 3.0f;  dy0 = y_iter - 1.0f;  dz0 = z_iter - 2.0f;

            float d1 = sqrtf(powf(dx1, 2) + powf(dy1, 2) + powf(dz1, 2));
            float d0 = sqrtf(powf(dx0, 2) + powf(dy0, 2) + powf(dz0, 2));
            float predicted_iter = d1 - d0;                         // predicted measurements in each iteration based on X_state
            float error_iter = measurement - predicted_iter;        // innovation term based on X_state

            // debug
            // error_iter = 0.1f;
            float e_y = error_iter;
            if ((d0 != 0.0f) && (d1 != 0.0f)){
                h[KC_STATE_X] = (dx1 / d1 - dx0 / d0);  
                h[KC_STATE_Y] = (dy1 / d1 - dy0 / d0); 
                h[KC_STATE_Z] = (dz1 / d1 - dz0 / d0);
                // e_y = linalg.inv(Q_c).dot(err_uwb)
                if (fabsf(Q_chol - 0.0f) < 0.0001f){
                    e_y = error_iter / 0.0001f;
                }
                else{ 
                    e_y = error_iter / Q_chol;
                }
                // e_x = inv(Ppr_c) * (error_x), here error_x = x_err
                // Problem: after deon mat_inv, Pc matrix becomes eye(9) !!!
                // Reason: arm_mat_inverse_f32() overwrites the source matrix !!!
                // https://community.arm.com/developer/tools-software/tools/f/keil-forum/32946/cmsis-dsp-matrix-inverse-problem
                matrixcopy(KC_STATE_DIM, KC_STATE_DIM, tmp1, P_chol);
                // in order to keep P_chol
                mat_inv(&tmp1m, &Pc_inv_m);                          // Pc_inv_m = inv(Pc_m) = inv(P_chol)
                mat_mult(&Pc_inv_m, &x_errm, &e_x_m);                  // e_x_m = Pc_inv_m.dot(x_errm) 

                // compute w_x, w_y --> weighting matrix
                // Since w_x is diagnal matrix, directly compute the inverse
                for (int state_k = 0; state_k < KC_STATE_DIM; state_k++){
                    GM_state(e_x[state_k], &wx_inv[state_k][state_k]);
                    wx_inv[state_k][state_k] = (float)1.0 / wx_inv[state_k][state_k];
                }

                // get w_y
                float w_y=0.0;
                GM_UWB(e_y, &w_y);

                // rescale covariance matrix P and Q
                mat_mult(&Pc_m, &wx_invm, &Pc_w_invm);       // Pc_w_invm = P_c.dot(linalg.inv(w_x))

                // rescale P matrix
                mat_mult(&Pc_w_invm, &Pc_tran_m, &P_w_m);        // P_w_m = Pc_w_invm.dot(Pc_tran_m) = P_c.dot(linalg.inv(w_x)).dot(P_c.T)
                // debug
                // matrixcopy(STATE_DIM, STATE_DIM, P_w, P);
                // rescale Q matrix
                float Q_w = 0.0f;
                if (fabsf(w_y - 0.0f) < 0.0001f){
                    Q_w = (Q_chol * Q_chol) / 0.0001f;
                }
                else{
                    Q_w = (Q_chol * Q_chol) / w_y;
                }
                // ====== INNOVATION COVARIANCE ====== //
                // debug
                // Q_w = tdoa->stdDev * tdoa->stdDev;
                // matrixcopy(STATE_DIM, STATE_DIM, P_w, P); 
                // H is a row vector
                mat_trans(&H, &HTm);
                mat_mult(&P_w_m, &HTm, &PHTm);     // PHTm = P_w.dot(H.T). The P is the updated P_w 

                float HPHR = Q_w;                  // HPH' + Q.            The Q is the updated Q_w 
                for (int i=0; i<KC_STATE_DIM; i++) {  // Add the element of HPH' to the above
                    HPHR += h[i]*PHTd[i];          // this obviously only works if the update is scalar (as in this function)
                }

                // ====== MEASUREMENT UPDATE ======
                // Calculate the Kalman gain and perform the state update
                for (int i=0; i<KC_STATE_DIM; i++) {
                    Kw[i] = PHTd[i]/HPHR;               // rescaled kalman gain = (PH' (HPH' + R )^-1) with the updated P_w and Q_w
                    //[Note]: The 'error' here is the innovation term based on x_check
                    x_err[i] = Kw[i] * error_check;           // error state for next iteration
                    X_state[i] = this->S[i] + x_err[i];  // convert to nominal state
                }
                // update P_iter matrix and Q matrix for next iteration
                matrixcopy(KC_STATE_DIM, KC_STATE_DIM, P_iter, P_w);
                Q_iter = Q_w;
            }
        }

        // ---------------------------------- After n iterations --------------------------------------- //
        // P = P_iter =P_w, arm matrix: Pm = P_w_m
        // Q = Q_iter = Q_w

        // Call the kalman update function with P,K and R=Q_iter
        kalmanCoreUpdateWithPKR(this, &HTm, &Kwm, error_check, Q_iter);

    } 
  }
  tdoaCount++;
}