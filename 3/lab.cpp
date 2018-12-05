#include <vector>
#include <cmath>
#include <iomanip>
#include <fstream>

#include "methods.h"
#include "constants.h"
#include "vec.h"
#include "funcs.h"

vec f(double t, const vec& x) {
	const double q_0 = q_nasos;
	const double p_0 = P_atm;
	double p_1 = x[0];
	double q_i = x[1];
	double p_i = x[2];
	double v   = x[3];
	double pos = x[4]; // x 
	double p_j = x[5];
	double q_j = x[6];
	double p_2 = x[7];
	double q_1 = x[8];

	if (v != 0) {
		v = v;
	}

	bool is_swap = pos < x_swap;

	vec result(9);

	if (is_swap)
		result[0] = Phi(q_0-q_i, p_1);
	else
		result[0] = Phi(q_1-q_i, p_1);

	result[1] = G(p_1 - p_i - P_alpha(q_i), q_i);
	result[2] = (q_i-s_i*v)/K_elast;
	result[3] = 1/m*(p_i*s_i-p_j*s_j-nu*v);
	result[4] = v;
	result[5] = (s_j*v-q_j)/K_elast;
	result[6] = G(p_j - p_2 - P_alpha(q_j), q_j);
	result[7] = Phi(q_j-q_1, p_2);

	if (is_swap)
		result[8] = G(p_2 - p_0 - P_alpha(q_1), q_1);
	else
		result[8] = G(p_1 - p_0 - P_alpha(q_1), q_1);

	return result;
}

int main() {
	vec start(9);
	start[0] = P_atm;
	start[1] = 0;
	start[2] = P_atm;
	start[3] = 0;
	start[4] = 0;
	start[5] = P_atm;
	start[6] = 0;
	start[7] = P_atm;
	start[8] = 0;

	double a = 0;
	double b = 30;
	double n = 30000000;
	double h = (b - a) / n;
	auto result = solveDE_Runge_Kutta4<vec>(a, b, h, start, f);

	print_data("res.dat", result, 1000);	
}