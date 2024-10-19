# this file containc C_$i_X, C_$i_Y needed for plot, see docs there
load 'C_coords.gp'

set terminal pdf size 4,6
set output 'main-cap.pdf'

set multiplot layout 2,1


set xlabel 'y(0)'
set ylabel 'z(0)' # offset -2,0
unset colorbox

set xrange [-4:4]
set yrange [-0.001:0.001]

set ytics -0.001,0.0004,0.001 offset 0.25,-0.25
set xtics -4,1,4 offset -0.25,-0.25

set margin -2
set tmargin -10
set lmargin at screen 0.2

plot \
    'BOX_G.dat' u 1:3:2:4 w boxxy fs solid 0.50 fc rgb 'light-gray' notitle,\
    'PG-Pimages-srt.dat' u 1:3:2:4 w boxxy fs solid 0.50 fc rgb 'dark-gray' notitle,\
    'BOX_C_0.dat' u 1:3:2:4 w boxxy fs solid 0.50 fc rgb 'light-red' lc rgb 'red' notitle,\
    'BOX_C_1.dat' u 1:3:2:4 w boxxy fs solid 0.50 fc rgb 'web-green' lc rgb 'green' notitle,\
    'BOX_C_2.dat' u 1:3:2:4 w boxxy fs solid 0.50 fc rgb 'web-blue' lc rgb 'blue' notitle,\
    'PC_0-Pimages-srt.dat' u ($1+C_1_X):($3+C_1_Y):2:4 w boxxy fs solid 0.50 fc rgb 'light-red' notitle,\
    'PC_1-Pimages-srt.dat' u ($1+C_2_X):($3+C_2_Y):2:4 w boxxy fs solid 0.50 fc rgb 'web-green' notitle,\
    'PC_2-Pimages-srt.dat' u ($1+C_0_X):($3+C_0_Y):2:4 w boxxy fs solid 0.50 fc rgb 'web-blue' notitle,\

set xlabel 'close tail'
set ylabel 'far tail' offset -2.5,0
unset colorbox

set xrange [-1.2:1.2]
set yrange [-1.2:1.2]

set ytics -1,0.5,1 offset 0.25,-0.25
set xtics -1,0.5,1 offset -0.25,-0.25

set margin -2
set tmargin -10    
set lmargin at screen 0.2

plot \
    'BOX_G.dat' u 5:7:6:8 w boxxy ls 1 fs solid 0.5 fc rgb 'light-gray' notitle,\
    'PG-Pimages-srt.dat' u 5:7:6:8 w boxxy ls 1 fs solid 0.5 fc rgb 'dark-gray' notitle,\
    
