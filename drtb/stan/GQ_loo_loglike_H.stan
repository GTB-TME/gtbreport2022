/* the log-likelihood vector needed by loo */

/* --- new --- */

/* surveillance */
for( j in 1:N_L_new_rr){
  log_lik[k] =
    multinomial_lpmf( L_new_rr_N[j] |
                      (nprops[L_new_rr_t[j]+1][L_new_rr_id[j]]*m1234)');
  k = k+1;
 }

for(j in 1:N_L_new_x){
  log_lik[k] = multinomial_lpmf( L_new_x_N[j] |
                                 (nprops[L_new_x_t[j]+1][L_new_x_id[j]]*m1234)');
  k = k+1;
 }


/* data created in TD block */
for(j in 1:N_L_new_dh0){
  log_lik[k] =
    multinomial_lpmf( L_new_dh0_N_noH[j] |
                      (nprops[L_new_dh0_t[j]+1][L_new_dh0_id[j]]*m1234)');
  k = k+1;
 }

/* survey */
for(j in 1:N_Y_new_rr){
  log_lik[k] =
    normal_lpdf( Y_new_rr_M[j] |
                 log((nprops[Y_new_rr_t[j]+1][Y_new_rr_id[j]]*m1234)'),
           Y_new_rr_S[j]);
  k = k+1;
 }

for(j in 1:N_Y_new_x){
  log_lik[k] =
    normal_lpdf( Y_new_x_M[j] |
                 log((nprops[Y_new_x_t[j]+1][Y_new_x_id[j]]*m1234)'),
                 Y_new_x_S[j]);
  k = k+1;
 }

for(j in 1:N_Y_new_dh0){        /* 1&2, 3&4 */
  log_lik[k] =
    normal_lpdf( Y_new_dh0_M_noH[j] |
             log((nprops[Y_new_dh0_t[j]+1][Y_new_dh0_id[j]]*m1234)'),
                 Y_new_dh0_S_noH[j]);
  k = k+1;
 }

for(j in 1:N_Y_new_dh1h){       /* 3&4 */
  log_lik[k] =
    normal_lpdf( Y_new_dh1h_M_noH[j] |
             log((nprops[Y_new_dh1h_t[j]+1][Y_new_dh1h_id[j]]*m34k)'),
                 Y_new_dh1h_S_noH[j]);
  k = k+1;
 }

/* --- ret --- */
/* surveillance */
for( j in 1:N_L_ret_rr){
  log_lik[k] =multinomial_lpmf(L_ret_rr_N[j] |
                               (rprops[L_ret_rr_t[j]+1][L_ret_rr_id[j]]*m1234)');
  k = k+1;
 }

for(j in 1:N_L_ret_x){
  log_lik[k] = multinomial_lpmf(L_ret_x_N[j]|
                                (rprops[L_ret_x_t[j]+1][L_ret_x_id[j]]*m1234)');
  k = k+1;
 }


/* data created in TD block */
for(j in 1:N_L_ret_dh0){
  log_lik[k]=multinomial_lpmf(L_ret_dh0_N_noH[j]|
                              (rprops[L_ret_dh0_t[j]+1][L_ret_dh0_id[j]]*m1234)');
  k = k+1;
 }

/* survey */
for(j in 1:N_Y_ret_rr){
  log_lik[k] =
    normal_lpdf( Y_ret_rr_M[j] |
                 log((rprops[Y_ret_rr_t[j]+1][Y_ret_rr_id[j]]*m1234)'),
                 Y_ret_rr_S[j]);
  k = k+1;
 }

for(j in 1:N_Y_ret_x){
  log_lik[k] =
    normal_lpdf( Y_ret_x_M[j] |
                 log((rprops[Y_ret_x_t[j]+1][Y_ret_x_id[j]]*m1234)'),
                 Y_ret_x_S[j]);
  k = k+1;
 }


/* data created in TD block */
for(j in 1:N_Y_ret_dh0){        /* 1&2, 3&4 */
  log_lik[k] =
    normal_lpdf( Y_ret_dh0_M_noH[j] |
                 log((rprops[Y_ret_dh0_t[j]+1][Y_ret_dh0_id[j]]*m1234)'),
                 Y_ret_dh0_S_noH[j]);
  k = k+1;
 }

for(j in 1:N_Y_ret_dh1h){       /* 3&4 */
  log_lik[k] =
    normal_lpdf( Y_ret_dh1h_M_noH[j] |
             log((rprops[Y_ret_dh1h_t[j]+1][Y_ret_dh1h_id[j]]*m34k)'),
                 Y_ret_dh1h_S_noH[j]);
  k = k+1;
 }

