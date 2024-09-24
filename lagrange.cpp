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

    double a = 0.0;
    double b = 1.5;
    
    std::ifstream infile("lagrange_data.dat");
    if (!infile.good()) {
        std::cerr << "Ошибка: файл lagrange_data.dat не найден или пуст." << std::endl;
        return 1;
    }

    std::ofstream dataFile("lagrange_data.dat");
    for (int n = 1; n <= 15; n++) {
        double* x = grid_step(a, b, n); //аргументы генерируем


        for (double X = a; X <= b; X += (b - a) / n) {
            double Y = laGrange(x, X, n);
            dataFile << X << " " << Y << " " << n << "\n";
        }
        dataFile << "\n"; 
        delete[] x; 
    }
    dataFile.close();

    std::ofstream gnuplot("sd.gp");
    gnuplot << "set title 'График интерполяции Лагранжа'\n";
    gnuplot << "set xlabel 'x'\n";
    gnuplot << "set ylabel 'y'\n";
    gnuplot << "set grid\n";
    gnuplot << "set xrange [0:1.5]\n";
    gnuplot << "set yrange [-1.5:0.5]\n";

   /* gnuplot << "plot ";
    for (int i = 1; i <= 15; ++i) {
        if (i > 1) {
            gnuplot<< ", " << std::endl;
        }
        gnuplot << "'lagrange_data.dat' index " << i - 1 << " title 'lagrange(" << i << ")' with lines lw 1 lc " << i;
    }*/
    /*std::string colors[] = { "red", "blue", "green", "orange", "purple",
                            "cyan", "magenta", "brown", "pink",
                            "yellow", "black", "violet", "seagreen",
                            "gold", "navy", "darkgreen" };*/

    
    gnuplot << "plot ";
    for (int i = 1; i <= 15; ++i) {
        gnuplot << "'lagrange_data.dat' using 1:" << i + 1 << " with lines title 'lagrange(" << i << "), ";
    }
    

    gnuplot << "\n";
    gnuplot.close();
    system("gnuplot -p sd.gp");
    //delete[]x;
}
