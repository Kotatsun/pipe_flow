#!/bin/sh


gnuplot<<EOF
set term gif animate delay 0.5
set output "./result/u.gif"
set size square            # same side lengths for x and y
set xlabel 'i'             # x-axis
set ylabel 'j'             # y-axis
set xrange[0:10]           # i-grid min & max
set yrange[0:1]           # j-grid min & max
set size 1,1
set parametric
set palette defined (-0.5 'blue', 0 'red', 0.5 'yellow')
set nosurface              # do not show surface plot
unset ztics                # do not show z-tics
set pm3d at b              # draw with colored contour
set view 0,0               # view from the due north

#set title 'rho at ${f%.txt} step '
#list = system("ls flow*.txt")
do for[fn in system("ls flow*.txt")]{
    splot fn using 1:2:3 ti fn
}
#splot for[fn in system("ls flow*.txt")] fn using 1:2:4 ti fn
#splot '$f' using 1:2:4
set output
#set terminal wxt enhanced
quit
EOF