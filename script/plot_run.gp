pars="count ln_post ln_like fAccept type Fn I0 Fs log_q log_L phi0 r0 tE t_0 noise_lev invtemp" 
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

plot      word(labels,ic)."_v_0512k_4.dat"  ev 100 u 1:(i==9?qfn300($9):i==12?r0fn(column(i)):(i==13?10**$13:column(i))) w lp ti "2014-".word(labels,ic)."-v (s".ic.")",\
      word(labels,ic)."_b_0512k_4.dat"  ev 100 u 1:(i==9?qfn300($9):i==12?r0fn(column(i)):(i==13?10**$13:column(i))) w lp ti "2014-".word(labels,ic)."-b (s".ic.")",\
