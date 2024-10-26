#include "Header.h"
using namespace std;

/*
int main() {
	double at = 0.0;
	double bt = 2 * M_PI;
	ofstream dataFile, dataFile2;

	
	//dataFile2.open("trig_delta.txt");

	dataFile.open("func_fnt.txt");
	for (double n = 1000; n <= 1000; ++n) {
		vector <double> t(2 * n + 2);
		t[0] = 0;
		for (int i = 1; i <= 2 * n + 1; i++) {
			t[i] = (2 * M_PI * (i - 1)) / (2 * n + 1);
		}

		double A0 = 0;
		double sum = 0;
		for (int i = 1; i <= 2 * n + 1; i++) {
			sum += func1(t[i]);
		}
		A0 = sum / (2 * n + 1);

		vector <double> ak(n + 1);
		vector <double> bk(n + 1);

		for (double j = 1; j <= n; j++) {
			double suma = 0;
			double sumb = 0;
			for (int i = 1; i <= 2 * n + 1; i++) {
				suma += func1(t[i]) * cos(j * t[i]);
				sumb += func1(t[i]) * sin(j * t[i]);
			}
			ak[j] = 2 * suma / (2 * n + 1);
			bk[j] = 2 * sumb / (2 * n + 1);
		}

		int index = 0;
		double delta = 0;
		double maxdelta = 0;
		for (double j = 3.557; j <= 3.5589; j += 2*M_PI/ 1e5) { //до t[2*n+1]
			double fnt = A0;
			for (double i = 1; i <= n; i++) {
				fnt += ak[i] * cos(i * j) + bk[i] * sin(i * j);
			}
			//double gt = func1(j);
			delta = abs(fnt - func1(j));
			if (delta > maxdelta) maxdelta = delta;
			dataFile << j << " " << fnt << " " << func1(j) << "\n";
		}
		if (maxdelta <= 1e-6) {
			index = n;
			break;
		}
		//dataFile.close();
		cout << "delta[" << n << "] = " << maxdelta << "\n";
		//dataFile2<< maxdelta << "\n";
	}
	dataFile.close();
	//dataFile2.close();

	std::ofstream gnuplot("sd.gp");
	gnuplot << "set grid\n";
	gnuplot << "set xlabel 't'\n";

	//gnuplot << "set xrange [3:4]\n";
	//gnuplot << "set yrange [0:10]\n";
	gnuplot << "set ylabel 'F(t)'\n";
	gnuplot << "plot 'func_fnt.txt' using 1:2  with lines title 'fnt', 'func_fnt.txt' using 1:3 with lines title 'gt'  \n";

	//gnuplot << "set xlabel 'n'\n";
	//gnuplot << "set ylabel 'delta'\n";
	//gnuplot << "plot 'trig_delta.txt'  with lines title 'delta'\n";
	//gnuplot << "plot 'func_fnt.txt' using 1:2 with lines title 'fnt', 'func_fnt.txt' using 1:3 with lines title 'gt'  \n";
	gnuplot.close();
	system("gnuplot -p sd.gp");
}*/