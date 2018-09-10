#include "methods.h"

std::vector<double> solveDE_Euler_Explicit(double a, double b, double h, double ya, DE_F f) {
    int n = (b-a)/h;

    std::vector<double> result(n+1, 0);
    result[0] = ya;

    double t_n = a;
    for (int i = 1; i < n+1; ++i) {
        double& y_n1 = result[i];
        double& y_n = result[i-1];

        y_n1 = y_n + h*f(t_n, y_n);

		t_n += h;
    }

    result.erase(result.begin());

    return result;
}
