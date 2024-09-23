set title 'График интерполяции Лагранжа'
set xlabel 'x'
set ylabel 'y'
set grid
set xrange [0:1.5]
set yrange [-2:1]
set pointsize 2
plot for [i=0:14] 'lagrange_data.dat' index i title sprintf('lagrange(%d)', i) with lines lw 2 lc i+1
