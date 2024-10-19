# this file containc C_$i_X, C_$i_Y needed for plot, see docs there
load 'C_coords.gp'

set terminal pdf size 4,3
set output 'rossler-3period-main-coords.pdf'
set xlabel 'y(0)'
set ylabel 'z(0)' # offset -2,0
unset colorbox

set xrange [-4:4]
set yrange [-0.001:0.001]

set ytics -0.001,0.0004,0.001 offset 0.25,-0.25
set xtics -4,1,4 offset -0.25,-0.25

set margin -2
set tmargin -10

plot \
    'BOX_G.dat' u 1:3:2:4 w boxxy fs solid 0.50 fc rgb 'light-goldenrod' notitle,\
    'PG-Pimages-srt.dat' u 1:3:2:4 w boxxy fs solid 0.50 fc rgb 'goldenrod' notitle,\
    'BOX_C_0.dat' u 1:3:2:4 w boxxy fs solid 0.50 fc rgb 'light-red' notitle,\
    'BOX_C_1.dat' u 1:3:2:4 w boxxy fs solid 0.50 fc rgb 'light-green' notitle,\
    'BOX_C_2.dat' u 1:3:2:4 w boxxy fs solid 0.50 fc rgb 'light-blue' notitle,\
    'PC_0-Pimages-srt.dat' u ($1+C_1_X):($3+C_1_Y):2:4 w boxxy fs solid 0.50 fc rgb 'red' notitle,\
    'PC_1-Pimages-srt.dat' u ($1+C_2_X):($3+C_2_Y):2:4 w boxxy fs solid 0.50 fc rgb 'green' notitle,\
    'PC_2-Pimages-srt.dat' u ($1+C_0_X):($3+C_0_Y):2:4 w boxxy fs solid 0.50 fc rgb 'blue' notitle,\
    

