#include <iostream>
#include <vector>
#include <iomanip>
#define _USE_MATH_DEFINES
#include <math.h>
#include <fstream>
//double pi = cos(-1);

double func1(double x) {
    double buf1 = M_PI * x;
    double buf2 = cos(buf1);
    return x * x * buf2;
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

double laGrange(double* x, double X, int n)
{
    double L = 0.0;
    for (int i = 0; i <= n; i++)
    {
        double Pr = 1.0;
        for (int j = 0; j <= n; j++)
        {
            if (i != j)
            {
                Pr = Pr * (X - x[j]) / (x[i] - x[j]);
            }
        }

        L += func1(x[i]) * Pr;
    }
    return L;
}

double Newton(double* x, double X, int n)
{
    double y[16];
    for (int i = 0; i <= n; i++)
        y[i] = func1(x[i]);
    double Pr = 1;
    double newton = y[0];
    for (int j = 1; j <= n; j++)
    {
        for (int i = 0; i <= n - j; i++)
        {
            y[i] = (y[i + 1] - y[i]) / (x[i + j] - x[i]);
        }
        Pr *= (X - x[j - 1]);
        newton = Pr * y[0] + newton;
    }
    return newton;
}


int main()
{
    std::cout << "Hello World!\n";

    double a = 0;
    double b = 3/2;
    int n = 15;

    double* x = grid_step(a, b, n);//ёмкость n+1
    
    double* res1 = new double[n + 1]; 
    double* res2 = new double[n + 1];
    double* res3 = new double[n + 1];

    for (int i = 0; i <= n; ++i) {
        res1[i] = func1(x[i]);
    }
    for (int i = 0; i <= n; ++i) {
        res2[i] = laGrange(x, x[i], n);
    }
    for (int i = 0; i <= n; ++i) {
        res3[i] = Newton(x, x[i], n);
    }


    for (int i = 0; i <= n; ++i) {
        std::cout << std::setprecision(15) << res1[i] << " ";

    }
    std::cout << "\n";
    for (int i = 0; i <= n; ++i) {
        std::cout << std::setprecision(15) << res2[i] << " ";
    }
    std::cout << "\n";
    for (int i = 0; i <= n; ++i) {
        std::cout << std::setprecision(15) << res3[i] << " ";

    }
    std::cout << "\n";
    //std::cout << std::setprecision(9) << laGrange(x, x[2], n);
    std::cout << "\n";
    for (int i = 0; i <= n; ++i) {
        std::cout << func1(x[i]) - laGrange(x, x[i], n) << "\n";
    }
    std::cout << "\n";
    for (int i = 0; i <= n; ++i) {
        std::cout << func1(x[i]) - Newton(x, x[i], n) << "\n";
    }
    std::cout << "\n";
    for (int i = 0; i <= n; ++i) {
        std::cout << x[i] << " ";

    }
    std::cout << "\n";
    std::cout << cos(3/2)<<"\n";
    std::cout << func1(0.5);

    std::ofstream gnuplot("sd.gp");
    gnuplot << "plot cos(x)";
    gnuplot.close();
    system("gnuplot -p sd.gp");
    //delete[]x;
}
