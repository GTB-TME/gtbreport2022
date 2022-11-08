/* data likelihood, new cases, not considering H */
/* NOTE also needs different J */

/* surveillance */
for( j in 1:N_L_new_rr){
  L_new_rr_N[j] ~
    multinomial(nprops[L_new_rr_t[j]+1][L_new_rr_id[j]]');
 }

for(j in 1:N_L_new_x){
  L_new_x_N[j] ~
    multinomial(nprops[L_new_x_t[j]+1][L_new_x_id[j]]');
 }

/* data created in TD block */
for(j in 1:N_L_new_dh0){
  L_new_dh0_N_noH[j] ~
    multinomial(nprops[L_new_dh0_t[j]+1][L_new_dh0_id[j]]');
 }

/* survey */
for(j in 1:N_Y_new_rr){
  Y_new_rr_M[j] ~
    normal(log(nprops[Y_new_rr_t[j]+1][Y_new_rr_id[j],2]),
           Y_new_rr_S[j]);
 }

for(j in 1:N_Y_new_x){
  Y_new_x_M[j] ~
    normal(log(nprops[Y_new_x_t[j]+1][Y_new_x_id[j],2]),
           Y_new_x_S[j]);
 }

/* data created in TD block */
for(j in 1:N_Y_new_dh0){
  Y_new_dh0_M_noH[j] ~
    normal(log(nprops[Y_new_dh0_t[j]+1][Y_new_dh0_id[j],2]),
           Y_new_dh0_S_noH[j]);
 }

for(j in 1:N_Y_new_dh1h){
  Y_new_dh1h_M_noH[j] ~
    normal(log(nprops[Y_new_dh1h_t[j]+1][Y_new_dh1h_id[j],2]),
           Y_new_dh1h_S_noH[j]);
 }


