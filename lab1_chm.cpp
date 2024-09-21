#include <iostream>
#include <vector>
#include <iomanip>
#define _USE_MATH_DEFINES
#include <math.h>


double func1(double x) {
    return x * x * cos(M_PI * x);
}

double* grid_step(double a, double b, int n) {
    double h = (b - a) / n;
    //std::vector<double> x(n + 1);
    double* x = new double[n + 1]; //make_unique?

    x[0] = a;
    x[n] = b;

    for (int i = 0; i < n; ++i) {
        x[i+1] = x[i] + h;
    }
    /*for (int i = 0; i <= n; ++i) {
        std::cout << x[i]<< " ";
    }*/
    return x;

}

double laGrange(double* x, double X, int n)//x-массив будет огромный
{
    double L = 0;
    for (int i = 0; i <= n; i++)
    {
        double Pr = 1;
        for (int j = 0; j <= n; j++)
        {
            if (i != j)
            {
                Pr = Pr * (X - x[j]) / (x[i] - x[j]);
            }
        }

        L = L + func1(x[i]) * Pr;
    }
    return L;
}



int main()
{
    std::cout << "Hello World!\n";

    double a = 0;
    double b = 1.5;
    int n = 2;

    double* first = grid_step(a, b, n);//ёмкость n+1
    
    double* res1 = new double[n+1]; //n+1
    for (int i = 0; i <= n; ++i) {
        res1[i] = func1(first[i]);
    }

    for (int i = 0; i <= 2; ++i) {
        std::cout << res1[i] << " ";

    }
    std::cout << "\n";
    for (int i = 0; i <= 2; ++i) {
        std::cout << first[i] << " ";

    }
    //delete[]first;
}
