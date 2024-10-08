#include "Header.h"
 

double func1(double x) {
    return x * x * cos(M_PI * x);
}

double* grid_step(double a, double b, int n) {
    double h = (b - a) / n;
    double* x = new double[n + 1];

    x[0] = a;
    for (int i = 1; i <= n; ++i) {
        x[i] = x[i - 1] + h;
    }
    return x;
}


double laGrange(double* x, double X, int n)
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

double Newton(double* x, double X, int n)
{
    double y[16];
    for (int i = 0; i <= n; i++) {
        y[i] = func1(x[i]);
    }
    double Pr = 1;
    double newton = y[0];
    for (int j = 1; j <= n; j++)
    {
        for (int i = 0; i <= n - j; i++)
        {
            y[i] = (y[i + 1] - y[i]) / (x[i + j] - x[i]);//Формула связи
        }
        Pr *= (X - x[j - 1]);
        newton = Pr * y[0] + newton;
    }
    return newton;
}


/*
int main()
{
    std::cout << "Hello World!\n";
    std::ofstream dataFile;
    double a = 0.0;
    double b = 1.5;
    int N = 15;
    

    //1.1
    for (int n = 1; n <= N; n++) {
        dataFile.open("lagrange_data_" + std::to_string(n) + ".txt"); //интерполяиця лагранжа
        double* arg = grid_step(a, b, n); 

        for (double X = a; X <= b; X += (b - a) / 1e5) {
            double Y = laGrange(arg, X, n);
            dataFile << X << " " << Y << " " << n << "\n";
        }
        dataFile << "\n";
        delete[] arg;
        dataFile.close();
    }
     
    //1.2 1.3
    dataFile.open("lagrange_deltaN.txt"); //ошибка приближения для каждого н

    for (int n = 1; n <= N; n++)
    {
        double* x = grid_step(a, b, n);
        double delta = 0;
        double maxdelta = 0;
        for (double i = a; i <= b; i += (b-a)/(1e5))
        {
            delta = abs(func1(i) - laGrange(x, i, n));
            if (delta > maxdelta)
            {
                maxdelta = delta;
                
            }
        }
        //std::cout << "delta[" << n << "] = " << maxdelta << std::endl;
        dataFile << maxdelta << std::endl;
        delete[]x;
    }
    dataFile.close();
    int N0 = 10;

    
    //1.4 
    dataFile.open("laGrangeErrN0.txt");//ошибка приближения для оптимального N0 Лагранжем и исходной функции
    for (double i = a; i <= b; i += (b - a) / (1e5)) //тут стоит количество точек для интерполирования
    {
        dataFile << i << " " << abs(func1(i) - laGrange(x, i, N0)) << std::endl; //тут стоит оптимальная степень 
        //std::cout << "y[" << i << "] - L[" << i << "] = " << abs(func1(i) - laGrange(x, i, 15)) << std::endl;
    }
    dataFile.close();  


    //2.1
    double cx[15];
    double cy[15];

    dataFile.open("chebyshev.txt");
    for (int i = 0; i <= N0-1; i++)//i = 0..n0-1
    {
        cx[i] = ((b + a) / 2) + ((b - a) / 2) * (cos(M_PI * (2 * i + 1) / (2 * N0)));
        cy[i] = func1(cx[i]);
        dataFile << cx[i] << " " << cy[i] << std::endl;
    }
    dataFile.close();

    //2.2
    dataFile.open("chebyshev_lagrange.txt");//многочлен лагранжа степени н0 по узлам чебышева
    for (int i = 0; i < N; i++)
    {
        //std::cout << " x = " << cx[i] << " L = " << laGrange(cx, cx[i], 15) << std::endl;
        dataFile << cx[i] << " " << laGrange(cx, cx[i], N0) << std::endl;
    }
    dataFile.close();

    //2.3
    dataFile.open("chebyshev_deltaN.txt");//оценка погрешности приближения дельта н на неравномерной сетке
    double maxdelta = 0;
    double delta = 0;
    for (double i = a; i <= b; i += (b - a) / (1e5))
    {
        delta = abs(func1(i) - laGrange(cx, i, N0));
        if (delta > maxdelta)
        {
            maxdelta = delta;
        }
    }

    //std::cout << "delta[" << N0 << "] = " << maxdelta << std::endl;
    dataFile << maxdelta << std::endl;

    dataFile.close();

    double* x = grid_step(a, b, N0);

    //2.4
    dataFile.open("cheb_lagrange_grid.txt");//сравнение двух многочленов Лагранжа Ln0(x) на равномерной и неравномерной сетках
    for (double i = 0.0; i <= 1.5; i += (b - a) / (1e5))
    {
        dataFile << i << " " << laGrange(x, i, 10) << " " << i << " " << laGrange(cx, i, 10) << std::endl;
    }
    dataFile.close();

    //2.4.2
    dataFile.open("cheb242.txt");
    for (double i = a; i <= b; i += (b - a) / (1e5-15)) 
    {
        dataFile << i << " " << abs(func1(i) - laGrange(cx, i, N0)) << std::endl; 
    }
    dataFile.close();


    //3.1
    dataFile.open("newtonN0.txt");//ньютон для оптимального N0
    double *nx = grid_step(a, b, 10);
    double newt[11];
    for (double i = a; i <= b; i+=(b-a)/1e5)
    {
        double nt = Newton(nx, i, N0);
        double lt = laGrange(nx, i, N0);
        dataFile << i << " " << nt <<  " " << func1(i) << " " << lt << std::endl;
        //std::cout << "Newton[" << n << "] = " << newt[n] << std::endl;
    }
    dataFile.close();

    //3.2 
    dataFile.open("newtonErrN0.txt");//Ошибка приближения Ньютона для оптимального н
    for (double i = a; i <= b; i += (b - a) / (1e5 - N0))
    {
        dataFile << i << " " <<abs(func1(i) - Newton(x, i, N0)) << std::endl;
        //std::cout << "y[" << i << "] - N[" << i << "] = " << abs(func1(i) - Newton(x, i, 15)) << std::endl;
    }
    dataFile.close();


    std::ofstream gnuplot("sd.gp");
    gnuplot << "set grid\n";
    //1.1
    /*gnuplot << "set title 'График интерполяции Лагранжа'\n";
    gnuplot << "set xlabel 'x'\n";
    gnuplot << "set ylabel 'y'\n";
    gnuplot << "set xrange [0:1.5]\n";
    gnuplot << "set yrange [-1.25:0.35]\n";
    gnuplot << "set key outside\n";
    gnuplot << "set lmargin 1\n";
    gnuplot << "set rmargin 23\n";
    gnuplot << "set bmargin 2\n";
    gnuplot << "set tmargin 2\n";
    gnuplot << "plot ";
    for (int i = 1; i <= 15; ++i) {
        gnuplot << "'lagrange_data_" << i << ".txt' using 1:2 with lines title 'lagrange(" << i << ")',";
    }*/
    

    //1.2 1.3
    /*gnuplot << "set title 'ошибка приближения лагранжем'\n";
    gnuplot << "set ylabel 'n'\n";
    gnuplot << "set ylabel 'delta n'\n";
    //gnuplot << "set xrange [1:15]\n";
    gnuplot << "set yrange [-0.1:1.3]\n";
    gnuplot << "plot 'lagrange_deltaN.txt' with lp title 'delta n'\n";*/


    //вроде не нужно было делать 1.3
    /*gnuplot << "set title 'Лагранж от N0'\n"; 
    gnuplot << "set xlabel 'x'\n";
    gnuplot << "set ylabel 'LaGrange N0 (x)'\n";
    gnuplot << "set yrange [-1.25:0.25]\n";
    gnuplot << "plot 'laGrangeN0.txt' with lines title 'Lagrange Interpolation N0' lc rgb 'blue' \n";*/
    //gnuplot << "plot 'func1_dat.txt' using 1:2 with lines title 'y(x) = x*x*cos(pi*x)' lc rgb 'red' \n"; два графика в один сливаются


    //1.4
    /*gnuplot << "set title 'Ошибка приближения для оптимального N0 Лагранжем'\n";//разница на узлах с исх функцией = 0
    gnuplot << "set xlabel 'x'\n";
    gnuplot << "set ylabel 'Ошибка'\n";
    gnuplot << "plot 'laGrangeErrN0.txt' with lines title 'y(x) - L(x)'\n";*/
       

    //2.3 Не надо
    /*gnuplot << "set title 'Ошибка приближения Лагранжем на неравномерной сетке'\n";
    gnuplot << "set xlabel 'n'\n";
    gnuplot << "set ylabel 'delta n'\n";
    gnuplot << "plot 'chebyshev_deltaN.txt' with lines title 'delta n'\n";*/

    //2.4
    /*gnuplot << "set title 'РАзница значений Лагранжа на равномерной и неравномерной сетке'\n";
    gnuplot << "set xlabel 'x'\n";
    //gnuplot << "set xrange [0:1.5]\n";
    gnuplot << "set ylabel 'y'\n"; 
    gnuplot << "plot 'cheb_lagrange_grid.txt'  using 1:2 with lines title 'L n0(x)', 'cheb_lagrange_grid.txt' using 3:4 with lines title 'Ln0(cx)\n";*/

    /*gnuplot << "set title 'Ошибка приближения для оптимального N0 Лагранжем на неравномерной сетке'\n";//разница на узлах с исх функцией = 0
    gnuplot << "set xlabel 'x'\n";
    gnuplot << "set ylabel 'Ошибка'\n";
    gnuplot << "plot 'cheb242.txt' with lines title 'y(x) - L(cx)'\n"; */



    //3.1
    /*gnuplot << "set title 'Интерполяция Ньютона для N0'\n";
    gnuplot << "set xlabel 'x'\n";
    gnuplot << "set ylabel 'Newton(x)'\n";
    //gnuplot << "plot 'newtonN0.txt' with lines title 'Newton(x)'\n";
    gnuplot << "plot 'newtonN0.txt'  using 1:2 with lines title 'Newton(x)', 'newtonN0.txt' using 1:3 with lines title 'y(x)', \
        'newtonN0.txt' using 1:4 with lines title 'lagrange(x)'\n";*/


    //3.2
    /*gnuplot << "set title 'Ошибка приближения для оптимального N0 Ньютон'\n";
    gnuplot << "set xlabel 'x'\n";
    gnuplot << "set ylabel 'Ошибка'\n";
    //gnuplot << "set yrange [-0.5e-15:1.7e-15]\n";
    gnuplot << "plot 'newtonErrN0.txt' with lines title 'y(x) - Newton(x)'\n";*/


    //gnuplot.close();
    //system("gnuplot -p sd.gp");
    
//}
