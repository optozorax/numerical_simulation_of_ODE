#include <vector>
#include <cmath>
#include <iomanip>
#include <fstream>

#include "methods.h"
#include "constants.h"
#include "vec.h"
#include "funcs.h"

void p(double t, vec& x) {
	// Перегородка, останавливающая груз
	if (x[4] > x_stop_r) {
		x[4] = x_stop_r;
		x[3] = 0;
	}
}

vec f(double t, const vec& x) {
	static bool swapped = false;
	const double p_0 = P_atm;
	static double p_1, q_i, p_i, v, pos, p_j, q_j, p_2, q_1;
	static vec result(9);

	p_1 = x[0];
	q_i = x[1];
	p_i = x[2];
	v   = x[3];
	pos = x[4]; // x 
	p_j = x[5];
	q_j = x[6];
	p_2 = x[7];
	q_1 = x[8];

	if (pos > x_swap_r) swapped = true;
	if (pos < x_swap_l) swapped = false;

	if (!swapped) {
		double q_0 = q_nasos;
		result[0] = Phi(q_0-q_i, p_1);
		result[1] = G(p_1 - p_i - P_alpha(q_i), q_i);
		result[2] = (q_i-s_i*v)/K_elast;
		result[3] = 1/m*(p_i*s_i-p_j*s_j-nu*v);
		result[4] = v;
		result[5] = (s_j*v-q_j)/K_elast;
		result[6] = G(p_j - p_2 - P_alpha(q_j), q_j);
		result[7] = Phi(q_j-q_1, p_2);
		result[8] = G(p_2 - p_0 - P_alpha(q_1), q_1);
	} else {
		double q_0 = -q_nasos;
		result[0] = Phi(q_1 - q_i, p_1);
		result[1] = G(p_1 - p_i - P_alpha(q_i), q_i);
		result[2] = (q_i - s_i * v) / K_elast;
		result[3] = 1 / m * (p_i*s_i - p_j * s_j - nu * v);
		result[4] = v;
		result[5] = (s_j*v - q_j) / K_elast;
		result[6] = G(p_j - p_2 - P_alpha(q_j), q_j);
		result[7] = Phi(q_j - q_0, p_2);
		result[8] = G(p_0 - p_1 - P_alpha(q_1), q_1);
	}

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
	double h = 0.0000003;
	double n = (b - a) / h;
	int save_count = 3000;
	bool is_for_latex = true;

	auto result = solveDE_Runge_Kutta4<vec>(a, b, h, save_count, start, f, true);
	if (is_for_latex) {
		for (auto& i : result) {
			i.second.push_back(0);
			i.second.push_back(q_nasos);
			i.second[0] -= P_atm;
			i.second[2] -= P_atm;
			i.second[5] -= P_atm;
			i.second[7] -= P_atm;
		}
		print_data("res.dat", result, true, {"p_1", "q_i", "p_i", "v", "x", "p_j", "q_j", "p_2", "q_1", "p_0", "q_0"});
	}
	print_data("res_excels.dat", result);
}