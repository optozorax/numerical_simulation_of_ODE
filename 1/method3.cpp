#include "methods.h"

std::vector<double> solveDE_Runge_Kutta4(double a, double b, double h, double ya, DE_F f) {
	int n = (b-a)/h;
	return std::vector<double>(n, 0);
}