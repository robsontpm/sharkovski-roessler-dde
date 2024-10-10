#c3[0].mid in coordinates: {[-2.917645895512649, -2.917645895512648],[-2.875075463023885e-05, -2.875075463023712e-05]}
#c3[1].mid in coordinates: {[-0.1200557313775255, -0.1200557313775255],[0.0002211630245879565, 0.0002211630245879566]}
#c3[2].mid in coordinates: {[3.364824104343936, 3.364824104343937],[-2.896911885108674e-05, -2.896911885108457e-05]}

# those are the coordiates translatrions needed to move the resulting box
# to the respective center. This must be here, because of the way the program
# now returns the images in good coordinates (around those points)
# So, the sets are 0,0- centered on output.
# To put them in good perspective of the full set, we need to translate them.

C_0_X = -2.917645895512649
C_0_Y = -2.875075463023885e-05

C_1_X = -0.1200557313775255
C_1_Y = 0.0002211630245879565

C_2_X = 3.364824104343936
C_2_Y = -2.896911885108674e-05

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
    
