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


plot \
    'bin/plots/BOX_grid3.dat' u 1:3:2:4 w boxxy fs solid 0.50 fc rgb 'light-goldenrod' notitle,\
    'bin/plots/Pgrid3.dat' u 1:3:2:4 w boxxy fs solid 0.50 fc rgb 'goldenrod' notitle,\
    'bin/plots/BOX_c3_0.dat' u 1:3:2:4 w boxxy fs solid 0.50 fc rgb 'light-red' notitle,\
    'bin/plots/BOX_c3_1.dat' u 1:3:2:4 w boxxy fs solid 0.50 fc rgb 'light-green' notitle,\
    'bin/plots/BOX_c3_2.dat' u 1:3:2:4 w boxxy fs solid 0.50 fc rgb 'light-blue' notitle,\
    'bin/plots/Pc3_0.dat' u ($1+c3_1_X):($3+c3_1_Y):2:4 w boxxy fs solid 0.50 fc rgb 'red' notitle,\
    'bin/plots/Pc3_1.dat' u ($1+c3_2_X):($3+c3_2_Y):2:4 w boxxy fs solid 0.50 fc rgb 'green' notitle,\
    'bin/plots/Pc3_2.dat' u ($1+c3_0_X):($3+c3_0_Y):2:4 w boxxy fs solid 0.50 fc rgb 'blue' notitle,\

