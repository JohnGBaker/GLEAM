all: test cp-test

test:
	(export OMP_NUM_THREADS=4 ; time ../../gleam -seed=0.01203453 -nskip=20 -nevery=100 -nsteps=5000 -prop=4 -de_ni=100 -de_reduce_gamma=16 -pt=6 -pt_evolve_rate=0.02 -remap_r0 -log_tE -model_extra_noise -Fn_max=20.5 -gen_data=m4.dat -gen_data_col=2 -gen_data_err_col=4 m4-p6_test_1 ) > m4-p6_test_1.out
	@diff m4-p6_test_[01]_t0.dat ; if [ $$? -ne 0 ] ; then  echo $@ FAIL ; else echo $@ PASS ;fi

cp-test:
	(export OMP_NUM_THREADS=4 ; time ../../gleam -seed=0.01203453 -nskip=20 -nevery=100 -nsteps=5000 -prop=4 -de_ni=100 -de_reduce_gamma=16 -pt=6 -pt_evolve_rate=0.02 -remap_r0 -log_tE -model_extra_noise -Fn_max=20.5 -gen_data=m4.dat -gen_data_col=2 -gen_data_err_col=4 -checkp_at_step=2000 m4-p6_test_2 ) > m4-p6_test_2.p1.out
	(export OMP_NUM_THREADS=4 ; time ../../gleam -seed=0.01203453 -nskip=20 -nevery=100 -nsteps=5000 -prop=4 -de_ni=100 -de_reduce_gamma=16 -pt=6 -pt_evolve_rate=0.02 -remap_r0 -log_tE -model_extra_noise -Fn_max=20.5 -gen_data=m4.dat -gen_data_col=2 -gen_data_err_col=4 -checkp_at_step=4000  -restart_dir=step_2000-cp m4-p6_test_2 ) > m4-p6_test_2.p2.out
	(export OMP_NUM_THREADS=4 ; time ../../gleam -seed=0.01203453 -nskip=20 -nevery=100 -nsteps=5000 -prop=4 -de_ni=100 -de_reduce_gamma=16 -pt=6 -pt_evolve_rate=0.02 -remap_r0 -log_tE -model_extra_noise -Fn_max=20.5 -gen_data=m4.dat -gen_data_col=2 -gen_data_err_col=4 -restart_dir=step_4000-cp m4-p6_test_2 ) > m4-p6_test_2.p3.out
	@diff m4-p6_test_[02]_t0.dat ; if [ $$? -ne 0 ] ; then  echo $@ FAIL ; else echo $@ PASS ;fi


(setenv OMP_NUM_THREADS 4;time ../../src/gleam/gleam -checkp_at_step=2000 -seed=0.01203453 -nchains=1 -nskip=20 -nevery=100 -pt_swap_rate=0.10 -nsteps=5000 -prop=4 -de_ni=100 -de_reduce_gamma=16 -pt=64 -pt_evolve_rate=0.02 -pt_Tmax=1000000000 -remap_r0 -log_tE -tE_max=150 -Fn_max=20.5 -tcut=-600 -gen_data=mock-2014-0270_4_mock.dat -gen_data_col=2 -gen_data_err_col=4 m4-2014-0270_h_0005k_4 ) >& m4-2014-0270_h_0005k_4.out
(setenv OMP_NUM_THREADS 4;time ../../src/gleam/gleam -restart_dir=step_2000-cp -checkp_at_step=3333 -seed=0.01203453 -nchains=1 -nskip=20 -nevery=100 -pt_swap_rate=0.10 -nsteps=5000 -prop=4 -de_ni=100 -de_reduce_gamma=16 -pt -pt_evolve_rate=0.02 -pt_n=64 -pt_Tmax=1000000000 -remap_r0 -log_tE -tE_max=150 -additive_noise -Fn_max=20.5 -tcut=-600 -gen_data=mock-2014-0270_4_mock.dat -gen_data_col=2 -gen_data_err_col=4 m4-2014-0270_h_0005k_4 ) >& m4-2014-0270_h_0005k_4_p1.out
(setenv OMP_NUM_THREADS 4;time ../../src/gleam/gleam -restart_dir=step_3333-cp -seed=0.01203453 -nchains=1 -nskip=20 -nevery=100 -pt_swap_rate=0.10 -nsteps=5000 -prop=4 -de_ni=100 -de_reduce_gamma=16 -pt -pt_evolve_rate=0.02 -pt_n=64 -pt_Tmax=1000000000 -remap_r0 -log_tE -tE_max=150 -additive_noise -Fn_max=20.5 -tcut=-600 -gen_data=mock-2014-0270_4_mock.dat -gen_data_col=2 -gen_data_err_col=4 m4-2014-0270_h_0005k_4 ) >& m4-2014-0270_h_0005k_4_p2.out&
