set xtics auto

set macros
unset key
set view map
set xtics border in scale 0,0 mirror norotate  offset character 0, 0, 0 autojustify
set ytics border in scale 0,0 mirror norotate  offset character 0, 0, 0 autojustify
set ztics border in scale 0,0 nomirror norotate  offset character 0, 0, 0 autojustify
#set nocbtics
set rtics axis in scale 0,0 nomirror norotate  offset character 0, 0, 0 autojustify
set cblabel "magnification" 
set cbrange [*:10] noreverse nowriteback
set palette rgbformulae -21, -22, -23
plot  basename."_mmap.dat" using 1:2:($3) with image

