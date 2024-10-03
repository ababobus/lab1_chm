set grid
set title 'РАзница значений Лагранжа на равномерной и неравномерной сетке'
set xlabel 'x'
set xrange [0:1.5]
set ylabel 'y'
plot 'cheb_lagrange_grid.txt' with lines title '(LN0(cx)) - LN0(x))'
