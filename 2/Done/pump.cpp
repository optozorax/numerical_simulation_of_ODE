#include "pump.h"

pump::pump() {

	X.resize(n + 1, 0);
	V.resize(n + 1, 0);
	Pi.resize(n + 1, 0);

	// Начальные значения
	X[0] = 0;
	V[0] = 0;
	Pi[0] = 1e5;
}

void pump::outputResults(char *fileName) {

	std::ofstream fout;
	fout.open(fileName);
	fout << "V\tPi\tX" << std::endl;
	fout << std::fixed << std::setprecision(5);
	for (int i = 0; i < n + 1; ++i) {
		fout << V[i] << "\t" << Pi[i] << "\t" << X[i] << std::endl;
	}
	fout.close();
}

void pump::solveDE_Euler_Explicit() {

	for (int i = 0; i < n; ++i) {

		V[i + 1] = V[i] + h*speed(Pi[i], V[i]);
		X[i + 1] = X[i] + h*xCoordinate(V[i]);
		Pi[i + 1] = Pi[i] + h*pressure(V[i]);
	}
}

void pump::solveDE_Euler_Nonexplicit() {

	for (int i = 0; i < n; ++i) {

		V[i + 1] = V[i] + h/2.0*(speed(Pi[i], V[i]) + speed(Pi[i] + h*speed(Pi[i], V[i]), V[i] + h*speed(Pi[i], V[i])));
		X[i + 1] = X[i] + h/2.0*(xCoordinate(V[i]) + xCoordinate(V[i] + h*speed(Pi[i], V[i])));
		Pi[i + 1] = Pi[i] + h/2.0*(pressure(V[i]) + pressure(V[i] + h*speed(Pi[i], V[i])));
	}
}

void pump::solveDE_Runge_Kutta4() {

	double k_1, k_2, k_3, k_4;
	for (int i = 0; i < n; ++i) {

		k_1 = speed(Pi[i], V[i]);
		k_2 = speed(Pi[i] + h/2.0*k_1, V[i] + h/2.0*k_1);
		k_3 = speed(Pi[i] + h/2.0*k_2, V[i] + h/2.0*k_2);
		k_4 = speed(Pi[i] + h*k_3, V[i] + h*k_3);
		V[i + 1] = V[i] + h/6.0*(k_1 + 2*k_2 + 2*k_3 + k_4);

		k_1 = xCoordinate(V[i]);
		k_2 = xCoordinate(V[i] + h/2.0*k_1);
		k_3 = xCoordinate(V[i] + h/2.0*k_2);
		k_4 = xCoordinate(V[i] + h*k_3);
		X[i + 1] = X[i] + h/6.0*(k_1 + 2*k_2 + 2*k_3 + k_4);

		k_1 = pressure(V[i]);
		k_2 = pressure(V[i] + h/2.0*k_1);
		k_3 = pressure(V[i] + h/2.0*k_2);
		k_4 = pressure(V[i] + h*k_3);
		Pi[i + 1] = Pi[i] + h/6.0*(k_1 + 2*k_2 + 2*k_3 + k_4);
	}
}