#include "Header.h"
using namespace std;
const double a = 0.0;
const double b = 1.5;

*/
void makeCoefficientsForPolynomial(std::vector<double>& coord_X, std::vector<double>& coord_Y, int N, std::vector<double>& a_i, std::vector<double>& b_i, std::vector<double>& c_i,
	std::vector<double>& d_i)
{
	for (int i = 1; i <= N - 1; ++i) {
		double h_i_prev = coord_X[i] - coord_X[i - 1];
		double h_i = coord_X[i] - coord_X[i - 1];
		b_i[i] = (coord_Y[i] - coord_Y[i - 1]) / h_i - 1. / 3 * h_i * (c_i[i + 1] + 2 * c_i[i]);
		d_i[i] = (c_i[i + 1] - c_i[i]) / (3 * h_i);
		a_i[i] = coord_Y[i - 1];
	}
	double h_n = coord_X[N] - coord_X[N - 1];
	b_i[N] = (coord_Y[N] - coord_Y[N - 1]) / h_n - 2. / 3 * h_n * c_i[N];
	d_i[N] = -c_i[N] / (3 * h_n);
	a_i[N] = coord_Y[N - 1];
}


void Progonka(const std::vector<double>& a_p, const std::vector<double>& b_p, const std::vector<double>& c_p, const std::vector<double>& d_p, std::vector<double>& c_i, int N)
{
	std::vector<double> dzeta(N + 2, 0);
	std::vector<double> eta(N + 2, 0);

	for (size_t i = 1; i <= N; ++i) {
		if (b_p[i] - a_p[i] * dzeta[i] != 0) {
			dzeta[i + 1] = c_p[i] / (b_p[i] - a_p[i] * dzeta[i]);
			eta[i + 1] = (a_p[i] * eta[i] - d_p[i]) / (b_p[i] - a_p[i] * dzeta[i]);
		}
	}
	for (size_t i = N; i >= 1; --i) {
		double value = dzeta[i + 1] * c_i[i + 1] + eta[i + 1];
		c_i[i] = dzeta[i + 1] * c_i[i + 1] + eta[i + 1];
	}
}

void makeCoefficientsForProgonka(std::vector<double>& a_p, std::vector<double>& b_p, std::vector<double>& c_p,
	std::vector<double>& d_p, std::vector<double>& coord_X, std::vector<double>& coord_Y, int N)
{
	for (int i = 2; i <= N; ++i) {
		double h_i_prev = coord_X[i - 1] - coord_X[i - 2];
		double h_i = coord_X[i] - coord_X[i - 1];
		a_p[i] = h_i_prev;
		b_p[i] = 2 * (h_i_prev + h_i);
		c_p[i] = h_i;
		d_p[i] = 3 * ((coord_Y[i] - coord_Y[i - 1]) /
			h_i - (coord_Y[i - 1] - coord_Y[i - 2]) / h_i_prev);
		//cout << a_p[i] << " " << b_p[i] << " " << c_p[i] << " " << d_p[i] << endl;
	}
}

/*
int main() {

	
	int N = 38;

	ofstream dataFile;
	double h = (b - a) / N;

	std::vector<double> coord_X, coord_Y;
	std::vector<double> a_i, b_i, c_i, d_i(N + 2, 0);
	std::vector<double> a_p, b_p, c_p, d_p(N + 2, 0);
	
	for (size_t i = 0; i <= N; i++)
	{
		coord_X.push_back(a + i * h);
		coord_Y.push_back(func1(a + i * h));
	}
	
	makeCoefficientsForProgonka(a_p, b_p, c_p, d_p, coord_X, coord_Y, N);
	Progonka(a_p, b_p, c_p, d_p, c_i, N);
	makeCoefficientsForPolynomial(coord_X, coord_Y, N, a_i, b_i, c_i, d_i);
	//polynomial_Y.push_back(a_i[1]);

	dataFile.open("spline_coordinates.txt");
	for (size_t i = 2; i <= N; ++i)
	{
		std::vector<double>polynomial_Y;
		

		double x_i_prev = a + (i - 1) * h;
		for (size_t j = 0; j <= 1e3; ++j)
		{
			double x_j = x_i_prev + (h / 1e3) * j;
			coord_X.push_back(x_j);
			polynomial_Y.push_back(a_i[i] + b_i[i] * (x_j - x_i_prev) + c_i[i] * pow((x_j - x_i_prev), 2) + d_i[i] * pow((x_j - x_i_prev), 3));
			coord_Y.push_back(func1(x_j));
			cout << coord_X[j] << " " << polynomial_Y[j] << endl;
			dataFile << coord_X[j] << " " << coord_Y[j] << endl;
		}
	}
	dataFile.close();
}*/
