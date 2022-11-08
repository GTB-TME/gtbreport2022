/* data likelihood for surveillance, new */

for(j in 1:N_Y_new_rr){
  Y_new_rr_M[j] ~
    normal(log((nprops[Y_new_rr_t[j]+1][Y_new_rr_id[j]]*m1234)'),
           Y_new_rr_S[j]);
 }

for(j in 1:N_Y_new_x){
  Y_new_x_M[j] ~
    normal(log((nprops[Y_new_x_t[j]+1][Y_new_x_id[j]]*m1234)'),
           Y_new_x_S[j]);
 }

for(j in 1:N_Y_new_dh0){
  Y_new_dh0_M[j] ~
    normal(log(nprops[Y_new_dh0_t[j]+1][Y_new_dh0_id[j],2:4]'),
           Y_new_dh0_S[j]);
 }

for(j in 1:N_Y_new_dh1h){
  Y_new_dh1h_M[j] ~
    normal(log((nprops[Y_new_dh1h_t[j]+1][Y_new_dh1h_id[j]]*m34k)'),
           Y_new_dh1h_S[j]);
 }

for(j in 1:N_Y_new_dh2){
  Y_new_dh2_M[j] ~
    normal(log((nprops[Y_new_dh2_t[j]+1][Y_new_dh2_id[j]]*m4k)'),
           Y_new_dh2_S[j]);
 }
