#define _USE_MATH_DEFINES

#include <vector>
#include <cmath>
#include <iomanip>
#include <fstream>
#include "methods.h"
#include "constants.h"
#include "vec.h"

//-----------------------------------------------------------------------------
double sign(double a) {
	if (a < 0)
		return -1;
	if (a > 0)
		return 1;
	if (a == 0)
		return 0;
}

//-----------------------------------------------------------------------------
double Phi(double qij, double p) {
	if (qij > 0 && p > 0 && p <= pAtm)
		return pow(p, 1.0 + 1.0 / gamma) * qij / C_cav;
	else
		return qij / C;
}

double P_alpha(double q) {
	if (mltRe*abs(q) < ReCnt)
		return r*q;
	else
		return r_kappa * pow(fabs(q), kappa) * sign(q);
}

double G(double dp, double q) {
	return B * sqrt(fabs(dp))*(F*pow(sqrt(fabs(dp / dzeta)), 3.0) - q);
}

//-----------------------------------------------------------------------------
vec f(double t, const vec& x) {
	const double q_1 = q_nasos;
	const double p_2 = P_atm;
	double p_1 = x[0];
	double q_2 = x[1];
	vec result(2);
	result[0] = Phi(q_1-q_2, p_1);
	result[1] = (p_1 - p_2) * fabs(q_1 - q_2); // Это работает в соответствии с ожиданиями
	//result[1] = G(p_1 - p_2 - P_alpha(q_1), q_1); // Это работает не так как ожидалось
	return result;
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

int main() {
	vec start(2);
	start[0] = 100000;
	start[1] = 0;
	double a = 0;
	double b = 0.1;
	double n = 1000;
	double h = (b - a) / n;
;	auto result = solveDE_Runge_Kutta4<vec>(a, b, h, start, f);

	std::ofstream fout("res.dat");
	//fout << "p_1\tq_2" << std::endl; // для LaTeX
	fout << std::setprecision(16);
	for (auto& i : result) {
		fout << i[0] << "\t" << i[1] << std::endl;
	}
	fout.close();
}