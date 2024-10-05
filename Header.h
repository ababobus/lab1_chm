#pragma once
#include <iostream>
#include <vector>
#include <iomanip>
#define _USE_MATH_DEFINES
#include <math.h>
#include <fstream>
#include <string>


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