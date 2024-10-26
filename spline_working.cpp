#include "Header.h"
using namespace std;
const double a = 0.0;
const double b = 1.5;


void Progonka(double* y, double* C, int n, double h) {
	double* alpha = new double[n + 1];
	double* beta = new double[n + 1];
	double* F = new double[n + 1];
	double z = 0;
	double c = 4 * h;
	double b = h;
	double a = h;
	int N = n - 1;
	for (int i = 1; i <= N; i++) {
		F[i] = (3 / h) * (y[i + 1] - 2 * y[i] + y[i - 1]);
	}
	alpha[1] = -c / b;
	beta[1] = F[1] / b;
	for (int i = 2; i <= N - 1; i++) {
		z = b + a * alpha[i - 1];
		alpha[i] = -c / z;
		beta[i] = (F[i] - a * beta[i - 1]) / z;
	}
	beta[N] = (F[N] - a * beta[N - 1]) / (b + a * alpha[N - 1]);
	C[n] = 0;
	C[N] = beta[N];
	for (int i = N - 1; i >= 1; i--) {
		C[i] = alpha[i] * C[i + 1] + beta[i];
	}
	C[0] = 0;
}
void coefPn(double* y, double* B, double* C, double* D, int n, double h) {
	for (int i = 1; i <= n - 1; i++) {
		D[i] = (C[i] - C[i - 1]) / (3 * h);
	}
	D[n] = -C[n - 1] / (3 * h);
	for (int i = 1; i <= n - 1; i++) {
		B[i] = ((y[i] - y[i - 1]) / h) - (h / 3) * (C[i] + 2 * C[i - 1]);//индексация 
	}
	B[n] = ((y[n] - y[n - 1]) / h) - (2 * h / 3) * C[n - 1];
}


double Pn(double xi, double* y, double* x, double* b, double* c, double* d, int n) {
	for (int i = 0; i <= n; i++) {
		if (xi >= x[i] && xi <= x[i + 1])
			return (y[i] + b[i + 1] * (xi - x[i]) + c[i] * pow((xi - x[i]), 2) + d[i + 1] * pow((xi - x[i]), 3));
	}
}



/*

int main()
{
	int n;
	std::ofstream dataFile;

	dataFile.open("spline.txt");
	for (int n = 1; n < 38; n++)
	{
		double h = (b - a) / n;
		double* x = new double[n + 1];
		double* y = new double[n + 1];
		for (int i = 0; i <= n; i++)
		{
			x[i] = i * h;
			y[i] = func1(x[i]);
		}
		double* B = new double[n + 1];
		double* C = new double[n + 1];
		double* D = new double[n + 1];
		Progonka(y, C, n, h);
		coefPn(y, B, C, D, n, h);
		double delta = 0.0;
		double maxDelta = 0.0;
		for (double i = 0; i <= b; i += h)
		{
			delta = abs(func1(i) - Pn(i, y, x, B, C, D, n));
			if (maxDelta <= delta)
				maxDelta = delta;
		}
		cout << "delta[" << n << "] = " << maxDelta;
		dataFile << n << " " << maxDelta << "\n";

		if (maxDelta <= 1e-6)
		{
			cout << " < 0.000001" << endl;
		}
	}
	dataFile.close();

	std::ofstream gnuplot("sd.gp");
	gnuplot << "set grid\n";
	gnuplot << "set title 'Ошибка приближения функции кубическим сплайном'\n";
	gnuplot << "set xlabel 'x'\n";
	gnuplot << "set ylabel 'Ошибка'\n";
	gnuplot << "plot 'spline.txt' with lines title 'y(x) - S(x)'\n";

	gnuplot.close();
	system("gnuplot -p sd.gp");
	return 0;
}*/
