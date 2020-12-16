#!/bin/sh


gnuplot<<EOF
set term gif animate delay 1.5
# set output "./result/u.gif"
set output "./result/the.gif"
set size square            # same side lengths for x and y
set xlabel 'x'             # x-axis
set ylabel 'y'             # y-axis
set xrange[0:20]           # i-grid min & max
set yrange[0:30]           # j-grid min & max
set size 1,1
set parametric
# set palette defined (310 'blue', 315 'cyan', 318 'green', 320 'yellow', 323 'red')
set palette defined(310"#0000ff",315"#0080ff",320"#00ffff",325"#00ff80",330"#80ff00",339"#ffff00",345"#ff8000",352"#ff0000")
set cbrange[310:352]
# set palette defined(310"#0000ff",315"#0080ff",320"#00ffff",325"#00ff80",330"#80ff00",335"#ffff00",340"#ff8000",345"#ff0000")
# set cbrange[310:345]
set nosurface              # do not show surface plot
unset ztics                # do not show z-tics
set pm3d at b              # draw with colored contour
set view 0,0               # view from the due north

# set title 'Temperature at ${f%.txt} step '
#list = system("ls flow*.txt")
do for[fn in system("ls output/flow*.txt")]{
    # splot fn using 1:2:3 ti fn
    splot fn using 1:2:4 ti fn
}
set output
quit
EOF