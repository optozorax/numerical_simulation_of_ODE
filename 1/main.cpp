#include <fstream>
#include "methods.h"
#include <iomanip>
#include <math.h>
const int _PRECISION_ = 5;

//-----------------------------------------------------------------------------
void writeTable(
	std::ofstream& fout,
	double a, double b,
	double h, double ya,
	DE_F f, DE_Solve solve
) {
	fout << "t      M1       M2       M3       Precise solve" << std::endl;
	auto result1 = solveDE_Euler_Explicit(a, b, h, ya, f);
	auto result2 = solveDE_Euler_Nonexplicit(a, b, h, ya, f);
	auto result3 = solveDE_Runge_Kutta4(a, b, h, ya, f);

	for (int i = 0; i < result1.size(); ++i) {
		a += h;
		fout << std::fixed << std::setprecision(3) << a << std::setprecision(_PRECISION_) << "  " << result1[i] << "  " << result2[i] << "  " << result3[i] << "  " << solve(a) << std::endl;
	}
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
double func(double t, double y) {
	return 2 * t * y;
}

//-----------------------------------------------------------------------------
double funcSolve(double t) {
	return exp(t*t);
}

//-----------------------------------------------------------------------------
int main() {
	std::ofstream fout;
	fout.open("table1.txt");
	writeTable(fout, 0, 1, 0.1, 1, func, funcSolve);
	fout.close();

	fout.open("table2.txt");
	writeTable(fout, 0, 1, 0.01, 1, func, funcSolve);
	fout.close();

	fout.open("table3.txt");
	writeTable(fout, 0, 1, 0.001, 1, func, funcSolve);
	fout.close();
}