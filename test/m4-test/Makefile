test:
	(export OMP_NUM_THREADS=4 ; time ../../gleam -seed=0.01203453 -nskip=20 -nevery=100 -nsteps=5000 -prop=4 -de_ni=100 -de_reduce_gamma=16 -pt=6 -pt_evolve_rate=0.02 -remap_r0 -log_tE -model_extra_noise -Fn_max=20.5 -gen_data=m4.dat -gen_data_col=2 -gen_data_err_col=4 m4-p6_test_1 ) > m4-p6_test_1.out
	@diff m4-p6_test_[01]_t0.dat ; if [ $$? -ne 0 ] ; then  echo FAIL ; else echo PASS ;fi

