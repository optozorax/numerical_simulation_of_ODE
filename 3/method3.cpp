#include "methods.h"

template<class T>
std::vector<T> solveDE_Runge_Kutta4(double a, double b, double h, T ya, std::function<T(double, T)> f) {
	int n = (b-a)/h;

	std::vector<T> result(n+1, 0);
	result[0] = ya;

	double t_n = a;
	for (int i = 1; i < n+1; ++i) {
		double t_n1 = t_n + h;
		T& y_n1 = result[i];
		T& y_n = result[i - 1];
		T k_n_1 = f(t_n, y_n);
		T k_n_2 = f(t_n + h / 2.0, y_n + h / 2.0 * k_n_1);
		T k_n_3 = f(t_n + h / 2.0, y_n + h /2.0 * k_n_2);
		T k_n_4 = f(t_n + h, y_n + h * k_n_3);

		y_n1 = y_n + h / 6.0 * (k_n_1 + 2.0 * k_n_2 + 2.0 * k_n_3 + k_n_4);

		t_n += h;
	}

	result.erase(result.begin());

	return result;
}