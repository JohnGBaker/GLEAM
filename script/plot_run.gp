pars="count ln_post ln_like fAccept type I0 Fs Fn log_q log_L r0 phi tE t_0 noise_lev" 
labels="July17_lc"
sinv(s,q0)=-1.+(q0+1.)/sqrt(1./s-1.)
qfn(s)=log10(sinv(s,1e5))
qfn1e3(s)=log10(sinv(s,1e3))
qfn300(s)=log10(sinv(s,300))
set ylabel word(pars,i)
r0fn(s)=1/sqrt(1/s-1)
tpk=2456700.
set pointsize 0.25
noise_mag(I0,noise_lev)=I0-2.5*log10(noise_lev)

plot      "gle_".word(labels,ic)."_v_0512k_4.dat"  ev 100 u 1:(i==9?qfn300($9):i==11?r0fn(column(i)):(i==8?1/0:(i==15?$8:(i==13?10**$13:column(i))))) w lp ti "2014-".word(labels,ic)."-v (s".ic.")",\
      "gle_".word(labels,ic)."_b_0512k_4.dat"  ev 100 u 1:(i==9?qfn300($9):i==11?r0fn(column(i)):(i==8?1/0:(i==15?$8:(i==13?10**$13:column(i))))) w lp ti "2014-".word(labels,ic)."-b (s".ic.")",\
      "gle_".word(labels,ic)."_c_0512k_4.dat"  ev 100 u 1:(i==9?qfn300($9):i==11?r0fn(column(i)):(i==8?1/0:(i==15?$8:(i==13?10**$13:column(i))))) w lp ti "2014-".word(labels,ic)."-c (s".ic.")",\
      "gle_".word(labels,ic)."_s_0512k_4.dat"  ev 100 u 1:(i==9?($9):i==11?r0fn(column(i)):(i==8?1/0:(i==15?$8:(i==13?10**$13:column(i))))) w lp ti "2014-".word(labels,ic)."-s (s".ic.")",\
      "gle_".word(labels,ic)."_r_0512k_4.dat"  ev 100 u 1:(i==9?($9):i==11?r0fn(column(i)):(i==8?1/0:(i==15?$8:(i==13?10**$13:column(i))))) w lp ti "2014-".word(labels,ic)."-r (s".ic.")",\
      "gle_".word(labels,ic)."_q_0512k_4.dat"  ev 100 u 1:(i==9?($9):i==11?r0fn(column(i)):(i==8?1/0:(i==15?$8:(i==13?10**$13:column(i))))) w lp ti "2014-".word(labels,ic)."-q (s".ic.")",\
      "gle_".word(labels,ic)."_p_0512k_4.dat"  ev 100 u 1:(i==9?($9):i==11?r0fn(column(i)):(i==8?1/0:(i==15?$8:(i==13?10**$13:column(i))))) w lp ti "2014-".word(labels,ic)."-p (s".ic.")",\
      "gle_".word(labels,ic)."_o_0512k_4.dat"  ev 100 u 1:(i==9?($9):i==11?r0fn(column(i)):(i==8?1/0:(i==15?$8:(i==13?10**$13:column(i))))) w lp ti "2014-".word(labels,ic)."-o (s".ic.")",\

plot\
      "gle_".word(labels,ic)."_p_0512k_4.dat"  ev 100 u 1:(i==9?($9):i==11?r0fn(column(i)):(i==8?1/0:(i==15?$8:(i==13?10**$13:column(i))))) w lp ti "2014-".word(labels,ic)."-p (s".ic.")",\
      "gle_".word(labels,ic)."_o_1024k_4.dat"  ev 100 u 1:(i==9?($9):i==11?r0fn(column(i)):(i==8?1/0:(i==15?$8:(i==13?10**$13:column(i))))) w lp ti "2014-".word(labels,ic)."-o+ (s".ic.")",\
      "gle_".word(labels,ic)."_n_1024k_4.dat"  ev 100 u 1:(i==9?($9):i==11?r0fn(column(i)):(i==8?1/0:(i==15?$8:(i==13?10**$13:column(i))))) w lp ti "2014-".word(labels,ic)."-n+ (s".ic.")",\
      "gle_".word(labels,ic)."_m_1024k_4.dat"  ev 100 u 1:(i==9?($9):i==11?r0fn(column(i)):(i==8?1/0:(i==15?$8:(i==13?10**$13:column(i))))) w lp ti "2014-".word(labels,ic)."-m+ (s".ic.")",\

