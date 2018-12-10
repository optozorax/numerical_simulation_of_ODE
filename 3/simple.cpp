#include <vector>
#include <cmath>
#include <iomanip>
#include <fstream>

#include "methods.h"
#include "constants.h"
#include "vec.h"
#include "funcs.h"

vec f(double t, const vec& x) {
	const double q_1 = q_nasos;
	const double p_2 = P_atm;
	double p_1 = x[0];
	double q_2 = x[1];
	vec result(2);
	result[0] = Phi(q_1-q_2, p_1);
	result[1] = G(p_1 - p_2 - P_alpha(q_2), q_2);
	return result;
}

int main() {
	vec start(2);
	start[0] = 100000;
	start[1] = 0;
	double a = 0;
	double b = 0.0002;
	double n = 1000000;
	double h = (b - a) / n;
	auto result = solveDE_Runge_Kutta4<vec>(a, b, h, start, f);

	//print_data("simple.dat", result, 1000);

	for (auto& i : result) i[0] -= P_atm;
	print_data("simple.dat", result, 1000, true, {"p", "q"}, a, b);
}