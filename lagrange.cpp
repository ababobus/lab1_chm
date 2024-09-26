#include "Header.h"
//double pi = cos(-1);

double func1(double x) {
    return x * x * cos(M_PI*x);
}

double* grid_step(double a, double b, int n) {
    double h = (b - a) / n;
    double* x = new double[n + 1]; 

    x[0] = a;
    //x[n] = b;

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

double Newton(double* x, double X, int n)
{
    double y[15+1];
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
    std::ofstream dataFile;
    double a = 0.0;
    double b = 1.5;
    double* x = grid_step(a, b, 15);

    
    for (int n = 1; n <= 15; n++) {
        dataFile.open("lagrange_data_" + std::to_string(n) + ".txt"); //интерполяиця лагранжа
        double* x = grid_step(a, b, n); 

        for (double X = a; X <= b; X += (b - a) / n) {
            double Y = laGrange(x, X, n);
            dataFile << X << " " << Y << " " << n << "\n";
        }
        dataFile << "\n";
        delete[] x;
        dataFile.close();
    }
     

    dataFile.open("lagrange_deltaN.txt"); //ошибка приближения для каждого н
    int N0 = 0; 
    for (int n = 1; n <= 15; n++)
    {
        double delta = 0;
        double maxdelta = 0;
        for (double i = a; i <= b; i += (b-a)/1e5)
        {
            delta = abs(func1(i) - laGrange(x, i, n));
            if (delta > maxdelta)
            {
                maxdelta = delta;
                N0 = n;
            }
        }
        //std::cout << "delta[" << n << "] = " << maxdelta << std::endl;
        dataFile << maxdelta << std::endl;
    }
    dataFile.close();

    dataFile.open("laGrangeN0.txt");//Лагранж для оптимального n 
    for (double i = a; i <= b; i += (b - a) / 15)
    {
        
        dataFile << laGrange(x, i, 15) << std::endl;
        //std::cout << "L[" << i << "] = " << laGrange(x, i, N0) << std::endl;
    }
    dataFile.close();



    dataFile.open("laGrangeErrN0.txt");//ошибка приближения для оптимального N0 Лагранжем и исходной функции
    for (double i = a; i <= b; i += (b - a) / 15)
    {
        dataFile << abs(func1(i) - laGrange(x, i, 15)) << std::endl;
        std::cout << "y[" << i << "] - L[" << i << "] = " << abs(func1(i) - laGrange(x, i, N0)) << std::endl;
    }
    dataFile.close();

    dataFile.open("newtonN0.txt");//ньютон для оптимального N0
    double* newt = new double[16];
    for (int n = 1; n <= 15; n++)
    {
        newt[n] = Newton(x, x[n], 15);
        dataFile << x[n] << " " << newt[n] << std::endl;
        //std::cout << "Newton[" << n << "] = " << newt[n] << std::endl;
    }
    dataFile.close();

    dataFile.open("newtonErrN0.txt");//Ошибка приближения Ньютона для оптимального н
    for (double i = a; i <= b; i += (b - a) / 15)
    {
        dataFile << abs(func1(i) - Newton(x, i, 15)) << std::endl;
        //std::cout << "y[" << i << "] - N[" << i << "] = " << abs(func1(i) - Newton(x, i, N0)) << std::endl;
    }
    dataFile.close();


    std::ofstream gnuplot("sd.gp");
    gnuplot << "set grid\n";
    /*gnuplot << "set title 'График интерполяции Лагранжа'\n";
    gnuplot << "set xlabel 'x'\n";
    gnuplot << "set ylabel 'y'\n";
    gnuplot << "set xrange [0:1.5]\n";
    gnuplot << "set yrange [-1.25:0.25]\n";
    gnuplot << "set key outside\n";
    gnuplot << "set lmargin 1\n";
    gnuplot << "set rmargin 23\n";
    gnuplot << "set bmargin 2\n";
    gnuplot << "set tmargin 2\n";
    gnuplot << "plot ";
    for (int i = 1; i <= 15; ++i) {
        gnuplot << "'lagrange_data_" << i << ".txt' using 1:2 with lines title 'lagrange(" << i << ")',";
    }*/
    
    
    /*gnuplot << "set title 'Ошибка приближения'\n";
    gnuplot << "set ylabel 'n'\n";
    gnuplot << "set ylabel 'delta n'\n";
    gnuplot << "set xrange [0:15]\n";
    gnuplot << "set yrange [-1:15]\n";
    gnuplot << "plot 'deltaN.txt' with lp title 'delta n'\n";*/

    /*gnuplot << "set title 'Лагранж от N0'\n";
    gnuplot << "set xlabel 'x'\n";
    gnuplot << "set ylabel 'LaGrange(x)'\n";
    gnuplot << "set xrange [0:1.5]\n";
    gnuplot << "set yrange [-1.25:0.25]\n";
    gnuplot << "plot 'laGrangeN0.txt' with lines title 'Lagrange Interpolation N0' lc rgb 'blue' \n";*/
    //gnuplot << "plot 'func1_dat.txt' using 1:2 with lines title 'y(x) = x*x*cos(pi*x)' lc rgb 'red' \n"; два графика в один сливаются

    gnuplot << "set title 'Ошибка приближения для оптимального N0 Лагранжем'\n";//разница на узлах с исх функцией = 0
    gnuplot << "set xlabel 'x'\n";
    gnuplot << "set ylabel 'Ошибка'\n";
    gnuplot << "plot 'laGrangeErrN0.txt' with lines title 'y(x) - L(x)'\n";
    
    /*gnuplot << "set title 'Интерполяция Ньютона для N0'\n";
    gnuplot << "set xlabel 'x'\n";
    gnuplot << "set ylabel 'Newton(x)'\n";
    gnuplot << "set xrange [0:1.5]\n";
    gnuplot << "set yrange [-1.25:0.25]\n";
    gnuplot << "plot 'newtonN0.txt' with lines title 'Newton(x)'\n";*/

    /*gnuplot << "set title 'Ошибка приближения для оптимального N0 Ньютон'\n";
    gnuplot << "set xlabel 'x'\n";
    gnuplot << "set ylabel 'Ошибка'\n";\
    gnuplot << "set yrange [-0.5e-15:1.7e-15]\n";
    gnuplot << "plot 'newtonErrN0.txt' with lines title 'y(x) - Newton(x)'\n";*/
    
    gnuplot.close();
    system("gnuplot -p sd.gp");
    
}
