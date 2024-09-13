#c3[0].mid in coordinates: {[-2.917645895512649, -2.917645895512648],[-2.875075463023885e-05, -2.875075463023712e-05]}
#c3[1].mid in coordinates: {[-0.1200557313775255, -0.1200557313775255],[0.0002211630245879565, 0.0002211630245879566]}
#c3[2].mid in coordinates: {[3.364824104343936, 3.364824104343937],[-2.896911885108674e-05, -2.896911885108457e-05]}

c3_0_X = -2.917645895512649
c3_0_Y = -2.875075463023885e-05

c3_1_X = -0.1200557313775255
c3_1_Y = 0.0002211630245879565

c3_2_X = 3.364824104343936
c3_2_Y = -2.896911885108674e-05


# bo musze przesunac do punktu referencyjnego, bo zapisuje coordynaty wzgledem nich a nie grid3...
# przepisałem je z danych wyżej (program gencoords je zwraca... Powinno byc automatycznie, ale na razie jest tak..)


set terminal pdf size 3,3
set output 'rossler-3period-C2-tail-coords.pdf'
set xlabel 'max Ai'
set ylabel 'max Xi' # offset -2,0
unset colorbox

set xrange [-1.2:1.2]
set yrange [-1.2:1.2]

set ytics -1,0.5,1 offset 0.25,-0.25
set xtics -1,0.5,1 offset -0.25,-0.25

set margin -2
set tmargin -10

#5:7:6:8

plot \
    'bin/plots/BOX_c3_0.dat' u 5:7:6:8 w boxxy fs solid 0.50 fc rgb 'light-red' notitle,\
    'bin/plots/Pc3_2.dat.csv' u 5:7:6:8 w boxxy fs solid 0.50 fc rgb 'blue' notitle,\
    
