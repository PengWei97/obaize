mpiexec -n 30  ~/projects/obaize/obaize-opt -i inp2_mech_only.i > 01.log & 
# --recover s3_t2_df_pf_mech_coupled/out_s3_t2_df_pf_mech_coupled_cp/LATEST
# tar -cvf - ex_s3_c2_pf_df_with_void1/* | pigz -9 -p 20 > ex_s3_c2_pf_df_with_void1.tgz