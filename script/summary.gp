if(exists("file")){
	filenamebase=file
} else {
    filenamebase = "ptmlc_2014-001_l_0256k_4"
}
chain=-1
pnames = "I0 Fs Fn log(q) log(L) log(r0) phi tE tmax"
npar=9
#pnames = "Ω_f  b  k  A_x  n_x  dt"
#set terminal postscript eps color enhanced font ",10"
multi=1
if(multi==1){
    set terminal pdf lw 4 enhanced font ",14" size 4*npar in,3*npar in
    set pointsize 3.0/npar
    set style line 2 lw 3.0/npar
    set style line 3 lw 3.0/npar
} else {
    set terminal pdf enhanced
    set pointsize 0.5
}
leadcol=5
#st0=58000
st0=200000
#st0=950000
#st0=500000
ievery=5
istart=3000
#st0=1250000
#ievery=100

sinv(s,q0)=-1.+(q0+1.)/sqrt(1./s-1.)
#qfn(s)=log10(sinv(s,1e5))
qfn(s)=s
#r0fn(s)=(1/sqrt(1/s-1))
r0fn(s)=log10(1/sqrt(1/s-1))

set grid
set key off
#ext=".out"
ext=".dat"

#name(ifile)="pmpc_xrun".ifile."_60k"
#name(ifile)="pmpc_xrun".ifile."_100k"
name(ifile)="pmpc_nrun".ifile."_100k"
#name(ifile)="pmpc_xrun".ifile."_5k"


    print "Making summary plot for ", filenamebase
    set output filenamebase.'_summary.pdf'
    if(multi==1){
	set multiplot layout npar,npar
    }
    do for [par2=1:npar] {
	print "plotting par ",par2
	do for [par1=1:npar] {
	    if ( par1 == par2 ) {
		set title word(pnames,par1).' posterior samples'
		#set key bottom right
		if( chain<0){
		    plot filenamebase.ext index 0 ev ievery::istart u ($1>st0?(par1==4?qfn(column(par1+leadcol)):par1==6?r0fn(column(par1+leadcol)):column(par1+leadcol)):1/0):(column(2)) ti "ch0",\
#		    "" index 1 ev ievery::istart u ($1>st0?(par1==4?qfn(column(par1+leadcol)):par1==6?r0fn(column(par1+leadcol)):column(par1+leadcol)):1/0):(column(2)) ti "ch1",\
		    "" index 2 ev ievery::istart u ($1>st0?(par1==4?qfn(column(par1+leadcol)):par1==6?r0fn(column(par1+leadcol)):column(par1+leadcol)):1/0):(column(2)) ti "ch2"
		} else {    
		  plot filenamebase.ext index chain ev ievery::istart u ($1>st0?(par1==4?qfn(column(par1+leadcol)):par1==6?r0fn(column(par1+leadcol)):column(par1+leadcol)):1/0):(column(2)) ti "ch".chain
		}
	    } else {	
		set title word(pnames,par1).' x '.word(pnames,par2)
		#set key bottom right
		if(chain<0){
			plot filenamebase.ext index 0 ev ievery::istart u ($1>st0?(par1==4?qfn(column(par1+leadcol)):par1==6?r0fn(column(par1+leadcol)):column(par1+leadcol)):1/0):((par2==4?qfn(column(par2+leadcol)):par2==6?r0fn(column(par2+leadcol)):column(par2+leadcol))) ti "ch0",\
#			"" index 1 ev ievery::istart u ($1>st0?(par1==4?qfn(column(par1+leadcol)):par1==6?r0fn(column(par1+leadcol)):column(par1+leadcol)):1/0):((par2==4?qfn(column(par2+leadcol)):par2==6?r0fn(column(par2+leadcol)):column(par2+leadcol))) ti "ch1",\
			"" index 2 ev ievery::istart u ($1>st0?(par1==4?qfn(column(par1+leadcol)):par1==6?r0fn(column(par1+leadcol)):column(par1+leadcol)):1/0):((par2==4?qfn(column(par2+leadcol)):par2==6?r0fn(column(par2+leadcol)):column(par2+leadcol))) ti "ch2"
		} else {
			plot filenamebase.ext index chain ev ievery*5::istart u ($1>st0?(par1==4?qfn(column(par1+leadcol)):par1==6?r0fn(column(par1+leadcol)):column(par1+leadcol)):1/0):((par2==4?qfn(column(par2+leadcol)):par2==6?r0fn(column(par2+leadcol)):column(par2+leadcol))) ti "ch".chain
		}
	    }
	}     
    }
    unset multiplot
#    set multiplot layout 2,2
#    print "plotting waves"
#    set title "1-sigma range"
#    if(chain<0){
#        print "chain = ",chain
#    	plot filenamebase."_1_sigma_samples_0.dat" u 1:2 w p ti "data" ls 1 ,"" u 1:3 w l ti "" ls 2, filenamebase."_1_sigma_samples_1.dat" u 1:3 w l ti "" ls 2, filenamebase."_1_sigma_samples_2.dat" u 1:3 w l ti "model" ls 2
#    	set title "best (each chain)"
#    	plot filenamebase."_best_0.dat" u 1:2 w p ti "data" ls 1,"" u 1:3 w l ti "" ls 3, filenamebase."_best_1.dat" u 1:3 w l ti "" ls 3, filenamebase."_best_2.dat" u 1:3 w l ti "best" ls 3, filenamebase."_best_0.dat" u 1:($3+$5) w l ti "" ls 2,"" u 1:($3-$5) w l ti "" ls 2, filenamebase."_best_1.dat" u 1:($3+$5) w l ti "" ls 2,"" u 1:($3-$5) w l ti "" ls 2, filenamebase."_best_2.dat" u 1:($3+$5) w l ti "400M" ls 2,"" u 1:($3-$5) w l ti "noise range" ls 2
#    } else {
#	plot filenamebase."_1_sigma_samples_".chain.".dat" u 1:2 w p ti "data" ls 1 ,"" u 1:3 w l ti "" ls 2
#    	set title "best (chain ".chain.")"
#    	plot filenamebase."_best_".chain.".dat" u 1:2 w p ti "data" ls 1,"" u 1:3 w l ti "" ls 3, "" u 1:($3+$5) w l ti "noise range" ls 2,"" u 1:($3-$5) w l ti "" ls 2
#    }       
#    #plot 0
    #plot 0
#    unset multiplot

set term X11 
print "done"


