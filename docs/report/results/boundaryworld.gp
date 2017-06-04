
projroot = '../../../'
files = projroot . 'build/analysis2/'
target = 'gen2/'

set key outside right box opaque spacing 1.2 height 1 width 1

set terminal epslatex input size 6in,3in

set style line 1 lc rgb '#dddddd' 
set style line 2 lc rgb '#2b83ba'
set style line 3 lc rgb '#fdae61'
set style line 4 lc rgb '#78ac70'
set style line 5 lc rgb '#d7191c'

set xlabel 'World Size $W$'
set ylabel 'Ghost Cell Fraction ${}^G/_{\left(C^{0.4562}\right)}$'

# data = "<( cat ".files."summary.dat | awk '{mean=0; for(i=2; i<14; i+=1) mean += $i / 12; print $1, mean }' )"

data = "<( cat ".target."summ2.dat | awk '{mean=0; for(i=2; i<14; i+=1) mean += $i / 12; print $1, mean }' )"


set output target."borderworldsize.tex"
# plot files.'summary.dat' u 1:2 w l ls 3 title 'Test Data', \
#     for [col=3:12] files.'summary.dat' u 1:col w l ls 3 notitle, \
#     data u 1:2 w l ls 2 lw 4 title 'Average'

plot target.'summ2.dat' u 1:2 w l ls 3 title 'Test Data', \
    for [col=3:12] target.'summ2.dat' u 1:col w l ls 3 notitle, \
    data u 1:2 w l ls 2 lw 4 title 'Average'


# unset output
# pause -1
