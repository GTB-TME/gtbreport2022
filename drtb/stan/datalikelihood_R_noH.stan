/* data likelihood, ret cases, not considering H */
/* NOTE also needs different J */

/* surveillance */
for( j in 1:N_L_ret_rr){
  L_ret_rr_N[j] ~
    multinomial(rprops[L_ret_rr_t[j]+1][L_ret_rr_id[j]]');
 }

for(j in 1:N_L_ret_x){
  L_ret_x_N[j] ~
    multinomial(rprops[L_ret_x_t[j]+1][L_ret_x_id[j]]');
 }

/* data created in TD block */
for(j in 1:N_L_ret_dh0){
  L_ret_dh0_N_noH[j] ~
    multinomial(rprops[L_ret_dh0_t[j]+1][L_ret_dh0_id[j]]');
 }

/* survey */
for(j in 1:N_Y_ret_rr){
  Y_ret_rr_M[j] ~
    normal(log(rprops[Y_ret_rr_t[j]+1][Y_ret_rr_id[j],2]),
           Y_ret_rr_S[j]);
 }

for(j in 1:N_Y_ret_x){
  Y_ret_x_M[j] ~
    normal(log(rprops[Y_ret_x_t[j]+1][Y_ret_x_id[j],2]),
           Y_ret_x_S[j]);
 }

/* data created in TD block */
for(j in 1:N_Y_ret_dh0){
  Y_ret_dh0_M_noH[j] ~
    normal(log(rprops[Y_ret_dh0_t[j]+1][Y_ret_dh0_id[j],2]),
           Y_ret_dh0_S_noH[j]);
 }

for(j in 1:N_Y_ret_dh1h){
  Y_ret_dh1h_M_noH[j] ~
    normal(log(rprops[Y_ret_dh1h_t[j]+1][Y_ret_dh1h_id[j],2]),
           Y_ret_dh1h_S_noH[j]);
 }
