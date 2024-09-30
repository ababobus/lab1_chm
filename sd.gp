set grid
set title 'Ошибка приближения для оптимального N0 Ньютон'
set xlabel 'x'
set ylabel 'Ошибка'
plot 'newtonErrN0.txt' with lines title 'y(x) - Newton(x)'
