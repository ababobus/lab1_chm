set grid
set title '������ ����������� ��� ������������ N0 ������'
set xlabel 'x'
set ylabel '������'
set yrange [-0.5e-15:1.7e-15]
plot 'newtonErrN0.txt' with lines title 'y(x) - Newton(x)'
