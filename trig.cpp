#include "Header.h"

const double A = 0;
const double B = 1.5;

double trig_interpol(double t, double* X, double* g, int n) {

	double A0 = 0;
	double Ak = 0;
	double Bk = 0;
	double sum = 0;
	double F = 0;
	for (int i = 0; i < 2 * n + 1; i++)
		sum += g[i];

	A0 = sum / (2 * n + 1);
	F += A0;

	for (int k = 1; k <= n; k++)
	{
		sum = 0;
		for (int i = 0; i < 2 * n + 1; i++)
			sum += g[i] * cos(k * X[i]);

		Ak = 2 * sum / (2 * n + 1);

		sum = 0;
		for (int i = 0; i < 2 * n + 1; i++)
			sum += g[i] * sin(k * X[i]);

		Bk = 2 * sum / (2 * n + 1);
		F += Ak * cos(k * t) + Bk * sin(k * t);
	}

	return F;

}


int main() {
	std::cout << "Hello world 2";
	std::ofstream dataFile;
	std::ofstream dataFile2;
	const int N = 10;
	double a = 0;
	double b = 2 * M_PI;
	double c = 1e3;

	dataFile.open("trig_delta.txt");
	dataFile2.open("trig_graphic.txt");
	for (int n = 1; n < 300; n++) {

		double h = (b - a) / n;

		double* t = new double[2 * n + 1];
		double* g = new double[2 * n + 1];

		for (int i = 0; i < 2 * n + 1; i++) {
			t[i] = 2 * M_PI * (i) / (2 * n + 1);
			g[i] = func1(t[i]);
		}

		double* t_rep = new double[c];
		double step = (b - a) / c;
		for (int i = 0; i < c - 1; i++)
			t_rep[i] = a + step * i;

		double MaxDelta = 0;
		for (int i = 0; i < c - 1; i++) {
			double delta = abs(trig_interpol(t_rep[i], t, g, n) - func1(t_rep[i]));//дельта 
			if (n == 100)
				dataFile2 << delta << std::endl;
			if (delta > MaxDelta)
				MaxDelta = delta;
		}
		std::cout << n % 10;
		if (MaxDelta < 1e-6) {
			std::cout << std::endl << "END\n" << "N = " << n << ", " << "  MaxDelta = " << MaxDelta << " < 2*10e3" << std::endl;
			break;
		}
		if (n % 10 == 0)
			std::cout << "\tN = " << n << ", " << "  MaxDelta = " << MaxDelta << std::endl;
		dataFile << MaxDelta << std::endl;
	}

	dataFile.close();
	dataFile2.close();


	std::ofstream gnuplot("sd.gp");

	gnuplot << "set grid\n";
	//gnuplot << "set xrange [0:2*pi]\n";
	gnuplot << "plot 'trig_graphic.txt' with lines title 'delta'\n";
	//gnuplot << "plot 'trig1.txt'  with lines title 'g(t)'\n";


	gnuplot.close();
	system("gnuplot -p sd.gp");

	return 0;

}