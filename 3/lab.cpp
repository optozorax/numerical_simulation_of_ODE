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
	const double q_0 = q_nasos;
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
		result[7] = Phi(q_0 - q_j, p_2);
		result[6] = G(p_2 - p_j - P_alpha(q_j), q_j);
		result[5] = -(s_j*(-v) - q_j) / K_elast;
		result[3] = -1/m*(-(p_i*s_i-p_j*s_j)-nu*(-v));
		result[4] = v;
		result[2] = -(q_i-s_i*(-v))/K_elast;
		result[1] = G(p_i - p_1 - P_alpha(q_i), q_i);
		result[0] = Phi(q_i - q_1, p_1);
		result[8] = G(p_1 - p_0 - P_alpha(q_1), q_1);
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
	int write_count = 3000;
	bool is_for_latex = true;

	std::ofstream fout("res.dat");
	if (is_for_latex)
		fout << "t\tq_0\tp_0\tp_1\tq_i\tp_i\tv\tx\tp_j\tq_j\tp_2\tq_1" << std::endl;
	fout << std::setprecision(16);
	int count = 0;

	vec k_n_1, k_n_2, k_n_3, k_n_4, y_n = start, y_n1;
	double t_n = a, t_n1;
	for (int i = 1; i < n + 1; ++i) {
		t_n1 = t_n + h;
		k_n_1 = f(t_n, y_n);
		k_n_2 = f(t_n + h / 2.0, y_n + h / 2.0 * k_n_1);
		k_n_3 = f(t_n + h / 2.0, y_n + h / 2.0 * k_n_2);
		k_n_4 = f(t_n + h, y_n + h * k_n_3);

		y_n1 = y_n + h / 6.0 * (k_n_1 + 2.0 * k_n_2 + 2.0 * k_n_3 + k_n_4);

		p(t_n1, y_n1);

		t_n = t_n1;
		y_n = y_n1;

		if (i % (int(n) / write_count) == 0) {
			std::cout << "\r" << std::setprecision(2) << std::fixed << std::setw(10) << double(count * 100) / write_count << "%";
			if (is_for_latex) {
				auto a = y_n;
				a[0] -= P_atm;
				a[2] -= P_atm;
				a[5] -= P_atm;
				a[7] -= P_atm;
				fout << t_n << "\t" << q_nasos << "\t" << 0 << "\t" << a << std::endl;
			} else {
				fout << y_n << std::endl;
			}
			count++;
		}
	}

	fout.close();

	//double b = 50, n = 90000000;
	//double b = 50.0/90.0, n = 1000000;
	//double h = (b - a) / n;
	//auto result = solveDE_Runge_Kutta4<vec>(a, b, h, start, f, p);
	//print_data("res.dat", result, 1000);	
}