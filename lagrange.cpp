#include "Header.h"
//double pi = cos(-1);

double func1(double x) {
    double buf1 = M_PI * x;
    double buf2 = cos(buf1);
    return x * x * buf2;
}

double* grid_step(double a, double b, int n) {
    double h = (b - a) / n;
    double* x = new double[n + 1]; 

    x[0] = a;
    x[n] = b;

    for (int i = 0; i < n; ++i) {
        x[i + 1] = x[i] + h;
    }
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

/*double Newton(double* x, double X, int n)
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
}*/


int main()
{
    std::cout << "Hello World!\n";
    std::ofstream dataFile;
    double a = 0.0;
    double b = 1.5;

    dataFile.open("lagrange_data.txt"); //интерполяиця лагранжа
    for (int n = 1; n <= 15; n++) {
        double* x = grid_step(a, b, n); 

        for (double X = a; X <= b; X += (b - a) / n) {
            double Y = laGrange(x, X, n);
            dataFile << X << " " << Y << " " << n << "\n";
        }
        dataFile << "\n";
        delete[] x;
    }
    dataFile.close();

    dataFile.open("deltaN.txt"); //ошибка приближения
    int c = 100000;
    double h = (b - a) / c;
    int N0 = 0; 
    double* x = grid_step(a, b, 15);
    for (int n = 1; n <= 15; n++)
    {
        double delta = 0;
        double sdelta = 0;
        for (double i = a; i <= b; i += h)
        {
            delta = abs(func1(i) - laGrange(x, i, n));
            if (sdelta < delta)
            {
                sdelta = delta;
                N0 = n;
            }
        }
        //std::cout << "Delta[" << n << "] = " << sdelta << std::endl;
        dataFile << sdelta << std::endl;
    }
    dataFile.close();

    dataFile.open("laGrangeN0.txt");//Лагранж для оптимального n 
    for (double i = a; i <= b; i += (b - a) / N0)
    {
        //std::cout << "L[" << i << "] = " << laGrange(x, i, N0) << std::endl;
        dataFile << laGrange(x, i, N0) << std::endl;
    }
    dataFile.close();

    dataFile.open("laGrangeErrN0.txt");//ошибка приближения для оптимального N0 Лагранжем 
    for (double i = a; i <= b; i += (b - a) / N0)
    {
        dataFile << abs(func1(i) - laGrange(x, i, N0)) << std::endl;//погрешность 
        std::cout << abs(func1(i) - laGrange(x, i, N0)) << std::endl;
    }
    dataFile.close();




    std::ofstream gnuplot("sd.gp");
    gnuplot << "set grid\n";
    //gnuplot << "set title 'График интерполяции Лагранжа'\n";
    //gnuplot << "set xlabel 'x'\n";
    //gnuplot << "set ylabel 'y'\n";
    //gnuplot << "set xrange [0:1.5]\n";
    //gnuplot << "set yrange [-1.5:0.5]\n";
    /*gnuplot << "plot ";
    for (int i = 0; i < 15; ++i) {
        gnuplot << "'lagrange_data.dat' using 1:" << i + 2 << " with lines title 'lagrange(" << i << "), ";
    }*/
    
    /*gnuplot << "set title 'Ошибка приближения'\n";
    gnuplot << "set ylabel 'n'\n";
    gnuplot << "set ylabel 'Delta'\n";
    gnuplot << "set xrange [0:15]\n";
    gnuplot << "set yrange [-1:15]\n";
    gnuplot << "plot 'deltaN.txt' with lp title 'Delta'\n";*/

    /*gnuplot << "set title 'Лагранж от n оптимального'\n";
    gnuplot << "set xlabel 'x'\n";
    gnuplot << "set ylabel 'L(x)'\n";
    gnuplot << "set xrange [0:15]\n";
    gnuplot << "plot 'laGrangeN0.txt' with lines title 'Lagrange Interpolation'\n";*/

    /*gnuplot << "set title 'Ошибка приближения для оптимального N0 Лагранжем'\n";//по нулям
    gnuplot << "set xlabel 'x'\n";
    gnuplot << "set ylabel 'Ошибка'\n";
    gnuplot << "plot 'laGrangeErrN0.txt' with lines title 'Ошибка приближения'\n";*/
    
    gnuplot.close();
    system("gnuplot -p sd.gp");
    
}
