set terminal png size 800,600
set output 'xy.png'
plot 'history.txt' using 1:2:4:5 with boxxyerrorbars
set output 'xz.png'
plot 'history.txt' using 1:3:4:6 with boxxyerrorbars
set output 'yz.png'
plot 'history.txt' using 2:3:5:6 with boxxyerrorbars
