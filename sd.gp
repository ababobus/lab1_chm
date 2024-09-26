set grid
set title 'График интерполяции Лагранжа'
set xlabel 'x'
set ylabel 'y'
set xrange [0:1.5]
set yrange [-1.25:0.25]
set key outside
set lmargin 1
set rmargin 23
set bmargin 2
set tmargin 2
plot 'lagrange_data_1.txt' using 1:2 with lines title 'lagrange(1)','lagrange_data_2.txt' using 1:2 with lines title 'lagrange(2)','lagrange_data_3.txt' using 1:2 with lines title 'lagrange(3)','lagrange_data_4.txt' using 1:2 with lines title 'lagrange(4)','lagrange_data_5.txt' using 1:2 with lines title 'lagrange(5)','lagrange_data_6.txt' using 1:2 with lines title 'lagrange(6)','lagrange_data_7.txt' using 1:2 with lines title 'lagrange(7)','lagrange_data_8.txt' using 1:2 with lines title 'lagrange(8)','lagrange_data_9.txt' using 1:2 with lines title 'lagrange(9)','lagrange_data_10.txt' using 1:2 with lines title 'lagrange(10)','lagrange_data_11.txt' using 1:2 with lines title 'lagrange(11)','lagrange_data_12.txt' using 1:2 with lines title 'lagrange(12)','lagrange_data_13.txt' using 1:2 with lines title 'lagrange(13)','lagrange_data_14.txt' using 1:2 with lines title 'lagrange(14)','lagrange_data_15.txt' using 1:2 with lines title 'lagrange(15)',