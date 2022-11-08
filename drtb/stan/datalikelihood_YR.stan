/* data likelihood for surveillance, ret */

for(j in 1:N_Y_ret_rr){
  Y_ret_rr_M[j] ~
    normal(log((rprops[Y_ret_rr_t[j]+1][Y_ret_rr_id[j]]*m1234)'),
           Y_ret_rr_S[j]);
 }

for(j in 1:N_Y_ret_x){
  Y_ret_x_M[j] ~
    normal(log((rprops[Y_ret_x_t[j]+1][Y_ret_x_id[j]]*m1234)'),
           Y_ret_x_S[j]);
 }

for(j in 1:N_Y_ret_dh0){
  Y_ret_dh0_M[j] ~
    normal(log(rprops[Y_ret_dh0_t[j]+1][Y_ret_dh0_id[j],2:4]'),
           Y_ret_dh0_S[j]);
 }

for(j in 1:N_Y_ret_dh1r){
  Y_ret_dh1r_M[j] ~
    normal(log((rprops[Y_ret_dh1r_t[j]+1][Y_ret_dh1r_id[j]]*m24k)'),
           Y_ret_dh1r_S[j]);
 }

for(j in 1:N_Y_ret_dh1h){
  Y_ret_dh1h_M[j] ~
    normal(log((rprops[Y_ret_dh1h_t[j]+1][Y_ret_dh1h_id[j]]*m34k)'),
           Y_ret_dh1h_S[j]);
 }

for(j in 1:N_Y_ret_dh2){
  Y_ret_dh2_M[j] ~
    normal(log((rprops[Y_ret_dh2_t[j]+1][Y_ret_dh2_id[j]]*m4k)'),
           Y_ret_dh2_S[j]);
 }
