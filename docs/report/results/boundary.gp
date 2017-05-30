
projroot = '../../../'
files = projroot . 'build/analysis/'
target = 'gen/'

set key outside right box opaque spacing 1.2 height 1 width 4

set terminal epslatex input size 6in,3in

# unset key

set style line 1 lc rgb '#dddddd' 
set style line 2 lc rgb '#2b83ba' lw 2
set style line 3 lc rgb '#fdae61' lw 2
set style line 4 lc rgb '#78ac70' lw 2
set style line 5 lc rgb '#d7191c' lw 2
set style line 6 lc rgb '#b324bd' lw 2

set xlabel 'Cell Count $\log_2{C}$'
set ylabel 'Ghost Count $\log_2{G}$'

# set label 1 at  2.55, 0.002598 '$W=2$';
# set label 2 at  2.55, 0.006885 '$W=4$';
# set label 3 at  2.55, 0.013261 '$W=8$';
# set label 4 at  2.55, 0.022926 '$W=16$';

set size ratio -1

set output target."mbrot-scale.tex"
plot \
    target.'mbrot-total.dat' u (log($1)/log(2)):(log($23)/log(2)) w l ls 6 title '$W=32$', \
    for [col=4:1:-1] target.'mbrot-total.dat' \
        u (log($1)/log(2)):(log(column(2**col))/log(2)) w l ls (col+1) title '$W='.(2**col).'$'

# pause -1
