# this file containc C_$i_X, C_$i_Y needed for plot, see docs there
load 'C_coords.gp'

set terminal pdf size 4,3
set output 'rossler-3period-X0-tail-coords.pdf'
set xlabel 'close tail'
set ylabel 'far tail' # offset -2,0
unset colorbox

set xrange [-1.2:1.2]
set yrange [-1.2:1.2]

set ytics -1,0.5,1 offset 0.25,-0.25
set xtics -1,0.5,1 offset -0.25,-0.25

set margin -2
set tmargin -10

plot \
    'BOX_G.dat' u 5:7:6:8 w boxxy ls 1 fs solid 0.5 fc rgb 'light-goldenrod' notitle,\
    'PG-Pimages-srt.dat' u 5:7:6:8 w boxxy ls 1 fs solid 0.5 fc rgb 'goldenrod' notitle,\
    
