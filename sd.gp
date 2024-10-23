set grid
set xlabel 't'
set xrange [3:4]
set yrange [0:10]
set ylabel 'F(t)'
plot 'func_fnt.txt' using 1:2  with lines title 'fnt', 'func_fnt.txt' using 1:3 with lines title 'gt'  
