
projroot = '../../../'
files = projroot . 'build/report/ghosts/'
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


set palette maxcolors 3
set palette model RGB defined (\
    0 '#CCEBFF', \
    1 '#FDD7B1', \
    2 '#C8FFBF')
unset colorbox

set output target."ghosts.tex"
plot files.'background.dat' with image pixels notitle, \
        files.'mesh.init.dat'w l ls 1 notitle, \
        files.'gh.r0.b1.dat'w l ls 2 title '$\mathbb{G}\vert_{p=0}$', \
        files.'gh.r0.b2.dat'w l ls 2 notitle, \
        files.'gh.r1.b0.dat'w l ls 3 title '$\mathbb{G}\vert_{p=1}$', \
        files.'gh.r1.b2.dat'w l ls 3 notitle, \
        files.'gh.r2.b0.dat'w l ls 4 title '$\mathbb{G}\vert_{p=2}$', \
        files.'gh.r2.b1.dat'w l ls 4 notitle

set output target."borders.tex"
plot files.'background.dat' with image pixels notitle, \
        files.'mesh.init.dat'w l ls 1 notitle, \
        files.'bd.r0.b1.dat'w l ls 2 title '$\mathbb{B}\vert_{p=0}$', \
        files.'bd.r0.b2.dat'w l ls 2 notitle, \
        files.'bd.r1.b0.dat'w l ls 3 title '$\mathbb{B}\vert_{p=1}$', \
        files.'bd.r1.b2.dat'w l ls 3 notitle, \
        files.'bd.r2.b0.dat'w l ls 4 title '$\mathbb{B}\vert_{p=2}$', \
        files.'bd.r2.b1.dat'w l ls 4 notitle

set output target."borderline.tex"
plot files.'background.dat' with image pixels notitle, \
        files.'mesh.init.dat'w l ls 1 notitle, \
        files.'bd.r0.dat'w l ls 2 title '$\mathbb{B}\vert_{p=0}$', \
        files.'bd.r1.dat'w l ls 3 title '$\mathbb{B}\vert_{p=1}$', \
        files.'bd.r2.dat'w l ls 4 title '$\mathbb{B}\vert_{p=2}$'

set output target."ghosts-r0.tex"
plot files.'background.dat' with image pixels notitle, \
        files.'mesh.init.dat'w l ls 1 notitle, \
        files.'bd.r0.dat'w l ls 2 title '$\mathbb{B}$', \
        files.'gh.r0.b1.dat'w l ls 3 lw 8 title '$\mathbb{G}_1$', \
        files.'gh.r0.b2.dat'w l ls 4 lw 8 title '$\mathbb{G}_2$', \
    
unset output
