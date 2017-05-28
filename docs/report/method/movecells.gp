
projroot = '../../../'
files = projroot . 'build/report/movecells/'
target = ''

set terminal epslatex input size 6in,3in

set size ratio -1
# set xrange [0.04:0.4]
# set yrange [0.25:0.6]
set xrange [0.2:0.6]
set yrange [0.15:0.55]

# set xtics 0.1
# set ytics 0.1,0.1,0.6
unset xtics
unset ytics

# set key at 0.39,0.59 box opaque spacing 1.2 height 1
# set key left top at 0.61,0.54 box opaque spacing 1.2 height 1
set key outside right box opaque spacing 1.2 height 1


set style line 1 lc rgb '#dddddd' 
set style line 2 lc rgb '#2b83ba' lw 6
set style line 3 lc rgb '#fdae61' lw 6
set style line 4 lc rgb '#78ac70' lw 6
set style line 5 lc rgb '#d7191c' lw 6
set style line 12 lc rgb '#2b83da' lw 8
set style line 15 lc rgb '#ff191c' lw 8


set output target."init.tex"
plot files.'mesh.init.dat'w l ls 1 notitle, \
        files.'hilb.r0.init.dat'w l ls 2 title '$\mathbb{L}_0$', \
        files.'hilb.r1.init.dat'w l ls 3 title '$\mathbb{L}_1$', \
        files.'hilb.r2.init.dat'w l ls 4 title '$\mathbb{L}_2$'

# pause -1

set output target."overlap.tex"
plot files.'mesh.init.dat'w l ls 1 notitle, \
        files.'hilb.r0.init.dat'w l ls 2 title '$\mathbb{L}_0$', \
        files.'hilb.r1.init.dat'w l ls 3 title '$\mathbb{L}_1$', \
        files.'hilb.r2.init.dat'w l ls 4 title '$\mathbb{L}_2$', \
        files.'hilb.overlapl.dat'w l ls 5 lw 8 title 'Send from 0 to 1'

# pause -1

set output target."ghostnotify.tex"
plot files.'mesh.init.dat'w l ls 1 notitle, \
        files.'hilb.r0.init.dat'w l ls 2 title '$\mathbb{L}_0$', \
        files.'hilb.r1.init.dat'w l ls 3 title '$\mathbb{L}_1$', \
        files.'hilb.r2.init.dat'w l ls 4 title '$\mathbb{L}_2$', \
        files.'hilb.overlapl.dat'w l ls 12 title 'To Send', \
        files.'ghostnotify.dat'w l ls 5 title '$\vb{\hat N}_2$'

# pause -1

set output target."Tr.tex"
plot files.'mesh.init.dat'w l ls 1 notitle, \
        files.'hilb.r0.init.dat'w l ls 2 title '$\mathbb{L}_0$', \
        files.'hilb.r1.init.dat'w l ls 3 title '$\mathbb{L}_1$', \
        files.'hilb.r2.init.dat'w l ls 4 title '$\mathbb{L}_2$', \
        files.'hilb.overlapl.dat'w l ls 12 title 'To Send', \
        files.'tr.dat'w l ls 15 lw 5 title '$T_r$'

set output target."result.tex"
plot files.'mesh.init.dat'w l ls 1 notitle, \
        files.'hilb.r0.final.dat'w l ls 2 title '$\mathbb{L}_0$', \
        files.'hilb.r1.final.dat'w l ls 3 title '$\mathbb{L}_1$', \
        files.'hilb.r2.final.dat'w l ls 4 title '$\mathbb{L}_2$'

unset output
