set xtics auto
set term pdf
set out basename.".pdf"


set macros
trajname = basename."_d_traj.dat"
command = "head -n 1 ".trajname.""
print "command is ".command
print system(command)
titlevar = system(command)
print "title is ", titlevar
set title titlevar font ",9"
plot basename."_d_lcrv.dat" u 2:(-$3) ti "data" ,basename."_lcrv.dat" u 2:(-$3) ti "typical fit" lw 2 w l,"" u 2:($4-$3) ti "" lw 2 w l,"" u 2:(-$4-$3) ti "model 1-sigma" lt 3 lw 2 w l
#plot [-1500:300] "2014-0001_b_c1_164k_d_lcrv.dat" u 2:(-$3)ps 0.7 ,"2014-0001_b_c1_164k_lcrv.dat" u 2:(-$3) lw 2 w l,"" u 2:($4-$3) lw 2 w l,"" u 2:(-$4-$3) lt 3 lw 2 w l

set xr [-100:100];rep;set xr [*:*]

unset key
set view map
set xtics border in scale 0,0 mirror norotate  offset character 0, 0, 0 autojustify
set ytics border in scale 0,0 mirror norotate  offset character 0, 0, 0 autojustify
set ztics border in scale 0,0 nomirror norotate  offset character 0, 0, 0 autojustify
#set nocbtics
set rtics axis in scale 0,0 nomirror norotate  offset character 0, 0, 0 autojustify
#set xrange [ -0.500000 : 4.50000 ] noreverse nowriteback
#set yrange [ -0.500000 : 4.50000 ] noreverse nowriteback
#set xrange [-2 :2] noreverse nowriteback
#set yrange [-2 :2] noreverse nowriteback
#set xrange [2938430 :2938470] noreverse nowriteback
#set yrange [-20 :20] noreverse nowriteback
set cblabel "magnification" 
#set cbrange [ -5 : 5.00000 ] noreverse nowriteback
set cbrange [*:10] noreverse nowriteback
#set palette rgbformulae -7, 2, -7
set palette rgbformulae -21, -22, -23
set size square
#plot  basename."_mmap.dat" using 1:2:(log10($3)-1) with image
plot  basename."_mmap.dat" using 1:2:($3) with image, trajname u 3:4  lt 1 pointsize 0.1
#plot  basename."_mmap.dat" using 1:2:($3) with image, basename."_d_traj.dat" u 3:4  lt 0 pointsize 0.1
#replot basename."_traj.dat" u 3:4 lt 3 w l

