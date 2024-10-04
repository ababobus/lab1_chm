set grid
set title 'Ошибка приближения для оптимального N0 Лагранжем'
set xlabel 'x'
set ylabel 'Ошибка'
plot 'laGrangeErrN0.txt' with lines title 'y(x) - L(x)'
