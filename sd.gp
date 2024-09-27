set grid
set title 'Ошибка приближения для оптимального N0 Ньютон'
set xlabel 'x'
set ylabel 'Ошибка'
set yrange [-0.5e-15:1.7e-15]
plot 'newtonErrN0.txt' with lines title 'y(x) - Newton(x)'
