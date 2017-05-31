
projroot = '../../../'
files = projroot . 'build/pois/'
target = 'gen/'

# set key outside right box opaque spacing 1.2 height 1 width 4
unset key

set terminal epslatex input size 6in,3in

set style line 1 lc rgb '#dddddd' 
set style line 2 lc rgb '#2b83ba' lw 2
set style line 3 lc rgb '#fdae61' lw 2
set style line 4 lc rgb '#78ac70' lw 2
set style line 5 lc rgb '#d7191c' lw 2
set style line 6 lc rgb '#b324bd' lw 2

set palette maxcolors 100
set palette model RGB defined (\
    0 '#2b83ba', \
    1 '#fdae61', \
    2 '#d7191c')
unset colorbox

# set xlabel 'Cell Count $\log_2{C}$'
# set ylabel 'Ghost Count $\log_2{G}$'

# set label 1 at  2.55, 0.002598 '$W=2$';
# set label 2 at  2.55, 0.006885 '$W=4$';
# set label 3 at  2.55, 0.013261 '$W=8$';
# set label 4 at  2.55, 0.022926 '$W=16$';

set size ratio -1
set ticslevel 0
set xtics 0.2
set ytics 0.2
set ztics 0.004,0.004,0.013
set zrange [0:0.012]
set view 55, 341, 1, 1
set dgrid3d 16,16
set pm3d interpolate 16,16


set output target."pois-splot.tex"
splot files.'P.final.dat' w pm3d ls 2

# pause -1
