#include "Header.h"

const double A = 0;
const double B = 1.5;

double trig_interpol() {
	for (int k = 1; k <= 15; ++k) {
		double* t = new double[2 * k + 2];
		t[0] = 0;
		t[2 * k + 2] = 2 * M_PI;
		
		double a0 = 0;

		double a[15];
		double b[15];

		double* x = grid_step(1, 2 * k + 1, 2 * k);
		for (int i = 1; i <= 2 * k + 1; i++) {
			t[i] = 2 * M_PI * (i - 1) / (2 * k + 1);

			a0 += func1(i);
			a[k] += func1(i) * cos(k * x[i]);
			b[k] += func1(i) * sin(k * x[i]);

		}
		a0 *= 1 / (2*k + 1);
		a[k] *= 2 / (2 * k + 1);
		b[k] *= 2 / (2 * k + 1);

	}

}


int main() {
	std::cout << "Hello world 2";
	std::ofstream dataFile;

	double t1 = 0.0;
	double t2 = 2.0 * M_PI;
	return 0;

}