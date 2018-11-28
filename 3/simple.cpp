#define _USE_MATH_DEFINES

#include <vector>
#include <cmath>
#include <iomanip>
#include <fstream>
#include "methods.h"

//-----------------------------------------------------------------------------
typedef std::vector<double> vec;
vec operator+(const vec& a, const vec& b) {
	#ifdef _DEBUG
	if (a.size() != b.size())
		throw std::exception();
	#endif
	vec result = a;
	for (int i = 0; i < b.size(); i++)
		result[i] += b[i];
	return result;
}

vec operator*(const vec& a, double b) {
	vec result = a;
	for (int i = 0; i < result.size(); i++)
		result[i] *= b;
	return result;
}

vec operator*(double b, const vec& a) {
	return operator*(a, b);
}

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
// Данные константы
const double q_nasos = 1.0e-3;	// расход воды слева
const double P_atm = 1e5;		// (атмосферное) давление справа
const double pi = 3.141592653589793238463;

// Константы, которые еще надо настроить
const double Si = 1;		// площадь цилиндра слева
const double Sj = 1;		// площадь цилиндра справа
const double w = 1;			// 
const double K_elast = 1;	// коэффициент упругости
const double m = 1;			// масса поршня
const double nu = 1;		// ню

const double d = 2e-3;		//
const double l = 5e-3;		//
const double v = 10;//M_PI * (d / 2)*(d / 2);
const double qij = 0;		//
const double gamma = 1.4;	//
const double pAtm = 1e5;	// атмосферное давление
const double ro = 1e3;		// плотность воды
const double ReCnt = 321;
const double Nu = 1e-6;		// число Нуссельта воды
const double dH = pi * d / (4 * Si);
const double kappa = 1;

// Функции, которые зависят только от констант

const double Es = 1260 * 1260 * 1000; //ro * pow(C, 2 * Si*n*d);	// модуль жёсткость воды ~10^9
const double C = v / Es;					// жёсткость воды
const double C_cav = pow(pAtm, 1 / gamma);	// жёсткость воды при кавитации
const double mltRe = dH / Nu * Si;
const double r = 12 * ro*Nu*l / (dH*dH*Si);
const double r_kappa = 0.1582 * ro * pow(Nu, 0.25)*l / (pow(dH, 1.25)*pow(Si, 1.75));

const double B = 1 / (l*sqrt(2 * ro));
const double F = Si * sqrt(2 * ro);
const double eps = 0.57 + 0.043*(1.1 - 1);
const double dzeta = 0.035*(pow(1 / eps - 1, 2));

//-----------------------------------------------------------------------------
double Phi(double qij, double p) {
	if (qij < 0)
		qij = qij;

	if (qij > 0 && p > 0 && p <= pAtm) {
		double _pow = 1.0 + 1.0 / gamma;
		p = pow(p, _pow);
		p = p * qij / C_cav;
		return p;
	} else
		return qij / C;
}

double P_alpha(double q) {
	if (mltRe*abs(q) < ReCnt)
		return r*q;
	else
		return r_kappa * pow(fabs(q), kappa) * sign(q);
}

double G(double dp, double q) {
	return B*sqrt(fabs(dp))*(F*pow(sqrt(fabs(dp / dzeta)), 3.0)-q) * sign(dp);
}

//-----------------------------------------------------------------------------
vec f(double t, const vec& x) {
	const double q_1 = q_nasos;
	const double p_2 = P_atm;
	double p_1 = x[0];
	double q_2 = x[1];
	vec result(2);
	result[0] = Phi(q_1-q_2, p_1);
	result[1] = G(p_1-p_2-P_alpha(q_2), q_2);
	if (result[1] != result[1]) {
		result[0] = Phi(q_1 - q_2, p_1);
		result[1] = G(p_1 - p_2 - P_alpha(q_2), q_2);
	}
	return result;
}

int main() {
	vec start(2);
	start[0] = 100000;
	start[1] = 0;
	auto result = solveDE_Euler_Explicit<vec>(0, 0.00001, 0.00001 /1000.0, start, f);

	std::ofstream fout("res.dat");
	//fout << "p_1\tq_2" << std::endl;
	fout << std::setprecision(16);
	for (auto& i : result) {
		//if (i[1] != i[1]) break;
		fout << i[0] << "\t" << i[1] << std::endl;
	}
	fout.close();
}