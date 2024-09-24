set title 'График интерполяции Лагранжа'
set xlabel 'x'
set ylabel 'y'
set grid
set xrange [0:1.5]
set yrange [-1.5:0.5]
plot 'lagrange_data.dat' using 1:2 with lines title 'lagrange(1), 'lagrange_data.dat' using 1:3 with lines title 'lagrange(2), 'lagrange_data.dat' using 1:4 with lines title 'lagrange(3), 'lagrange_data.dat' using 1:5 with lines title 'lagrange(4), 'lagrange_data.dat' using 1:6 with lines title 'lagrange(5), 'lagrange_data.dat' using 1:7 with lines title 'lagrange(6), 'lagrange_data.dat' using 1:8 with lines title 'lagrange(7), 'lagrange_data.dat' using 1:9 with lines title 'lagrange(8), 'lagrange_data.dat' using 1:10 with lines title 'lagrange(9), 'lagrange_data.dat' using 1:11 with lines title 'lagrange(10), 'lagrange_data.dat' using 1:12 with lines title 'lagrange(11), 'lagrange_data.dat' using 1:13 with lines title 'lagrange(12), 'lagrange_data.dat' using 1:14 with lines title 'lagrange(13), 'lagrange_data.dat' using 1:15 with lines title 'lagrange(14), 'lagrange_data.dat' using 1:16 with lines title 'lagrange(15), 
