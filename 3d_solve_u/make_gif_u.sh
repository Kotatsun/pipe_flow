#!/bin/sh


gnuplot<<EOF
set term gif animate delay 1.5
set output "./result/u.gif"
# set output "./result/the.gif"
set size square            # same side lengths for x and y
set xlabel 'x'             # x-axis
set ylabel 'y'             # y-axis
set xrange[0:40]           # i-grid min & max
set yrange[0:45]           # j-grid min & max
set size 1,1
set parametric
set palette defined (310 'blue', 310.5 'red', 311 'yellow')
# set cbrange[310:311]
set nosurface              # do not show surface plot
unset ztics                # do not show z-tics
set pm3d at b              # draw with colored contour
set view 0,0               # view from the due north

# set title 'Temperature at ${f%.txt} step '
#list = system("ls flow*.txt")
do for[fn in system("ls output/flow*.txt")]{
    splot fn using 1:2:3 ti fn
    # splot fn using 1:2:4 ti fn
}
set output
quit
EOF