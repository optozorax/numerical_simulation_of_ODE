#pragma once

#include <vector>
#include <string>

//-----------------------------------------------------------------------------
inline double sign(double a) {
	if (a < 0)
		return -1;
	if (a > 0)
		return 1;
	if (a == 0)
		return 0;
}

//-----------------------------------------------------------------------------
inline double Phi(double qij, double p) {
	if (qij < 0 && p > 0 && p <= pAtm)
		return pow(p, 1.0 + 1.0 / gamma) * qij / C_cav;
	else
		return qij / C;
}

inline double P_alpha(double q) {
	if (mltRe*abs(q) < ReCnt)
		return r*q;
	else
		return r_kappa * pow(fabs(q), kappa) * sign(q);
}

inline double G(double dp, double q) {
	return B * sqrt(fabs(dp))*(F*pow(sqrt(fabs(dp / dzeta)), 3.0)*sign(dp) - q);
}

//-----------------------------------------------------------------------------
template<class T>
void print_data(std::string file_name, const std::vector<T>& data, int write_count, bool is_for_latex = false, const std::vector<std::string>& latex_str = {}, double a = 0, double b = 0) {
	std::ofstream fout(file_name);
	if (is_for_latex) {
		fout << "t\t";
		for (auto& i : latex_str)
			fout << i << "\t";
		fout << std::endl;
	}
	fout << std::setprecision(16);
	double mul = (b - a) / write_count;
	int count = 0;
	for (auto& i : data) {
		if (count % (data.size()/write_count) == 0) {
			if (is_for_latex)
				fout << (a + count*mul) << "\t" << i << std::endl;
			else
				fout << i << std::endl;
		}
		count++;
	}
	fout.close();
}