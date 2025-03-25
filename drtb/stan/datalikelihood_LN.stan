/* data likelihood for surveillance, new */

for( j in 1:N_L_new_rr){
  L_new_rr_N[j] ~
    multinomial((nprops[L_new_rr_t[j]+1][L_new_rr_id[j]]*m1234)');
 }

for(j in 1:N_L_new_hr){
  L_new_hr_N[j] ~
    multinomial((nprops[L_new_hr_t[j]+1][L_new_hr_id[j]]*m1324)');
 }

for(j in 1:N_L_new_rm){
  /* NOTE new data order */
  L_new_rm_N[j] ~
    multinomial((nprops[L_new_rm_t[j]+1][L_new_rm_id[j]]*m4)');
 }

for(j in 1:N_L_new_x){
  L_new_x_N[j] ~
    multinomial((nprops[L_new_x_t[j]+1][L_new_x_id[j]]*m1234)');
 }

for(j in 1:N_L_new_dh0){
  L_new_dh0_N[j] ~
    multinomial(nprops[L_new_dh0_t[j]+1][L_new_dh0_id[j]]');
 }

for(j in 1:N_L_new_dh1){
  L_new_dh1_N[j] ~
    multinomial((nprops[L_new_dh1_t[j]+1][L_new_dh1_id[j]]*m24)');
 }


for(j in 1:N_L_new_dh2){
  L_new_dh2_N[j] ~
    multinomial((nprops[L_new_dh2_t[j]+1][L_new_dh2_id[j]]*m4)');
 }

