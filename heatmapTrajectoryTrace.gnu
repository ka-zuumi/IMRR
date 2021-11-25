set terminal pngcairo enhanced size 1200,600
set output "heatmapTrajectoryTrace.png"

inputfile="HBrCO2library/expcompareGradientsInGridtoNWChemGradients2_history3_20/data/interpolation.dat"
Nhistory=3

set multiplot
unset key
set border 3 lw 2
set xtics nomirror
set ytics nomirror

set tmargin at screen 0.95
set bmargin at screen 0.15

set grid x y lw 2
set lmargin at screen 0.05
set rmargin at screen 0.45
set xrange [1:11]
set yrange [2:12]

set object 1 rect from 6,6 to 7,7 fc rgb "light-green" fs border lc rgb "light-green"

set xlabel "r_{H-C} (A)"
set ylabel "r_{Br-O,min} (A)" offset char 1,0
plot inputfile u 1:2 w l lw 1 lc rgb "black",\
     inputfile u ($3>Nhistory?($1):1/0):2 w l lw 2 lc rgb "dark-green"


set tmargin at screen 0.75
set bmargin at screen 0.35

unset object
unset grid
set grid y lw 2
set lmargin at screen 0.55
set rmargin at screen 0.95
set xrange [0:10000]
set autoscale y

set xlabel "Time (Steps)"
set ylabel "Ninterpolation" offset char 1,0
plot inputfile u 0:3 w l lw 1 lc rgb "black",\
     inputfile u ($3>Nhistory?($0):1/0):3 w l lw 2 lc rgb "dark-green"
