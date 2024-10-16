#include "Header.h"
using namespace std;
const double ax = 0;
const double bx = 1.5;
double TrigonometricInterpol(double t, double* x, int n)
{
	double A0 = 0;
	double Ak = 0;
	double Bk = 0;
	double sum = 0;
	double result = 0;
	for (int i = 0; i < 2 * n + 1; i++)
		sum += func1(i);

	A0 = sum / (2 * n + 1);
	result += A0;

	for (int k = 1; k <= n; k++)
	{
		sum = 0;
		for (int i = 0; i < 2 * n + 1; i++)
			sum += func1(i) * cos(k * x[i]);

		Ak = 2 * sum / (2 * n + 1);

		sum = 0;
		for (int i = 0; i < 2 * n + 1; i++)
			sum += func1(i) * sin(k * x[i]);

		Bk = 2 * sum / (2 * n + 1);
		result += Ak * cos(k * t) + Bk * sin(k * t);
	}

	return result;
}
/*
int main() {
	double at = 0.0;
	double bt = 2 * M_PI;
	ofstream dataFile;


	dataFile.open("func_fnt.txt");
	for (double n = 1; n < 200; n++) {
		vector <double> t(2 * n + 2);
		
		for (int i = 1; i <= 2 * n + 1; i++) {
			t[i] = (2 * M_PI * (i-1)) / (2 * n + 1);//узлы 
		}

		double A0 = 0;
		double sum = 0;
		for (int i = 1; i <= 2 * n + 1; i++) {
			sum += func1(t[i]);
		}
		A0 = sum / (2 * n + 1);

		//vector <double> ak(n);
		//vector <double> bk(n);
		double* ak = new double[n+1];
		double* bk = new double[n+1];
		for (int j = 1; j <= n; j++) {
			double suma = 0;
			double sumb = 0;
			for (int i = 1; i <= 2 * n + 1; i++) {
				suma += func1(t[i]) * cos((j) * t[i]);
				sumb += func1(t[i]) * sin((j) * t[i]);
			}
			ak[j] = 2 / (2 * n + 1) * suma;
			bk[j] = 2 / (2 * n + 1) * sumb;
		}


		double fnt = A0;
		double delta = 0;
		double maxdelta = 0;
		for (double j = 0; j < 2 * M_PI; j += (2 * M_PI / 1e5)) {

			for (int k = 1; k <= n; k++) {

				fnt += ak[k] * cos(k * j) + bk[k] * sin(k * j);//ak[i - 1] * cos(i * t[j]) + bk[i - 1] * sin(i * t[j]);

				double gt = func1(j);
				delta = abs(fnt - gt);
				if (delta > maxdelta)
				{
					maxdelta = delta;

				}
				//if (maxdelta < 1e-6) ind = k;
				
			}

		}
		std::cout << "delta[" << n << "] = " << maxdelta << std::endl;
		dataFile << maxdelta << endl;
		
		
	}
	dataFile.close();

	dataFile.open("func_gt.txt");
	for (double j = 0; j < 2 * M_PI; j += 2 * M_PI / 1e5) {
		double gt = func1(j);
		dataFile << j << " " << gt << "\n";
	}
	dataFile.close();


	//dataFile.open("func_gt.txt");
	//for (int i = 1; i <= n; ++i) {
	//	double delta = 0;
	//	double maxdelta = 0;
	//	for (double j = 0; j < 2 * M_PI; j += 2 * M_PI / 1e5) {
	//		
	//		delta = abs(func1(j) - )
	//	}
	//}
	
	dataFile.close();

	std::ofstream gnuplot("sd.gp");
	gnuplot << "set grid\n";
	//gnuplot << "plot 'func_fnt.txt' with lines title 'fnt', ";
	//gnuplot << "'func_gt.txt' with lines title 'gt'\n";
	gnuplot << "plot 'func_fnt.txt' with lines title 'fnt'\n";
	gnuplot.close();
	system("gnuplot -p sd.gp");

	//double a = 0;
	//double b = 2 * M_PI;
	//ofstream dataFile;

	//dataFile.open("trig_delta.txt");
	//for (int n = 1; n <= 50; n++) {
	//	double* arg = new double[2 * n + 1];
	//	double* t = new double[2 * n + 1];
	//	double* g = new double[2 * n + 1];

	//	for (int i = 0; i < 2 * n + 1; i++) {
	//		t[i] = 2 * M_PI * (i) / (2 * n + 1);//узлы
	//		g[i] = func1(3 * i / (4 * M_PI)); //x=alpha*t+beta
	//		arg[i] = t[i] * 3 / (4 * M_PI);
	//		//dataFile << t[i] << " " << g[i] << endl;
	//	}

	//	double delta = 0;
	//	double maxdelta = 0;

	//	for (double i = 0; i < 2 * M_PI; i += (2 * M_PI) / (1000))
	//	{
	//		delta = abs(func1(3 * i / (4 * M_PI)) - TrigonometricInterpol(i, arg, n));
	//		if (delta > maxdelta)
	//		{
	//			maxdelta = delta;

	//			dataFile << i << " " << delta << endl;
	//		}
	//		//dataFile << maxdelta << std::endl;
	//		dataFile.close();

	//	}
	//	delete[] t;
	//	delete[] g;

	//	std::ofstream gnuplot("sd.gp");
	//	gnuplot << "set grid\n";
	//	gnuplot << "plot 'trig_delta.txt' with lines title 'delta'\n ";

	//	gnuplot.close();
	//	system("gnuplot -p sd.gp");
	//	return 0;
	//}
}*/
