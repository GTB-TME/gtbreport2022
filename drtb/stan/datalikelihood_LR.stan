/* data likelihood for surveillance, ret */

for( j in 1:N_L_ret_rr){
  L_ret_rr_N[j] ~
    multinomial((rprops[L_ret_rr_t[j]+1][L_ret_rr_id[j]]*m1234)');
 }

for(j in 1:N_L_ret_hr){
  L_ret_hr_N[j] ~
    multinomial((rprops[L_ret_hr_t[j]+1][L_ret_hr_id[j]]*m1324)');
 }

for(j in 1:N_L_ret_rm){
  /* NOTE data order */
  L_ret_rm_N[j] ~
    multinomial((rprops[L_ret_rm_t[j]+1][L_ret_rm_id[j]]*m4)');
 }

for(j in 1:N_L_ret_x){
  L_ret_x_N[j] ~
    multinomial((rprops[L_ret_x_t[j]+1][L_ret_x_id[j]]*m1234)');
 }

for(j in 1:N_L_ret_dh0){
  L_ret_dh0_N[j] ~
    multinomial(rprops[L_ret_dh0_t[j]+1][L_ret_dh0_id[j]]');
 }

for(j in 1:N_L_ret_dh1){
  L_ret_dh1_N[j] ~
    multinomial((rprops[L_ret_dh1_t[j]+1][L_ret_dh1_id[j]]*m24)');
 }


for(j in 1:N_L_ret_dh2){
  L_ret_dh2_N[j] ~
    multinomial((rprops[L_ret_dh2_t[j]+1][L_ret_dh2_id[j]]*m4)');
 }
