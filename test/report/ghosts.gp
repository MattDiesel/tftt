
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

set output "ghosts.tex"
plot 'background.dat' with image pixels notitle, \
        'mesh.init.dat'w l ls 1 notitle, \
        'gh.r0.b1.dat'w l ls 2 title '$\mathbb{G}\vert_{p=0}$', \
        'gh.r0.b2.dat'w l ls 2 notitle, \
        'gh.r1.b0.dat'w l ls 3 title '$\mathbb{G}\vert_{p=1}$', \
        'gh.r1.b2.dat'w l ls 3 notitle, \
        'gh.r2.b0.dat'w l ls 4 title '$\mathbb{G}\vert_{p=2}$', \
        'gh.r2.b1.dat'w l ls 4 notitle

set output "borders.tex"
plot 'background.dat' with image pixels notitle, \
        'mesh.init.dat'w l ls 1 notitle, \
        'bd.r0.b1.dat'w l ls 2 title '$\mathbb{B}\vert_{p=0}$', \
        'bd.r0.b2.dat'w l ls 2 notitle, \
        'bd.r1.b0.dat'w l ls 3 title '$\mathbb{B}\vert_{p=1}$', \
        'bd.r1.b2.dat'w l ls 3 notitle, \
        'bd.r2.b0.dat'w l ls 4 title '$\mathbb{B}\vert_{p=2}$', \
        'bd.r2.b1.dat'w l ls 4 notitle

set output "borderline.tex"
plot 'background.dat' with image pixels notitle, \
        'mesh.init.dat'w l ls 1 notitle, \
        'bd.r0.dat'w l ls 2 title '$\mathbb{B}\vert_{p=0}$', \
        'bd.r1.dat'w l ls 3 title '$\mathbb{B}\vert_{p=1}$', \
        'bd.r2.dat'w l ls 4 title '$\mathbb{B}\vert_{p=2}$'

set output "ghosts-r0.tex"
plot 'background.dat' with image pixels notitle, \
        'mesh.init.dat'w l ls 1 notitle, \
        'bd.r0.dat'w l ls 2 title '$\mathbb{B}$', \
        'gh.r0.b1.dat'w l ls 3 lw 8 title '$\mathbb{G}_1$', \
        'gh.r0.b2.dat'w l ls 4 lw 8 title '$\mathbb{G}_2$', \
    
unset output
