set grid
set title '»нтерпол€ци€ Ќьютона дл€ N0'
set xlabel 'x'
set ylabel 'Newton(x)'
plot 'newtonN0.txt'  using 1:2 with lines title 'Newton(x)', 'newtonN0.txt' using 1:3 with lines title 'y(x)',         'newtonN0.txt' using 1:4 with lines title 'lagrange(x)'
