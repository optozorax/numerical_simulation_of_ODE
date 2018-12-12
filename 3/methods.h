#pragma once

#include <vector>
#include <functional>

//-----------------------------------------------------------------------------
/** Явный метод Эйлера для решения обыкновенных ДУ: y'=f(t, y). */
template<class T>
std::vector<std::pair<double, T>> solveDE_Euler_Explicit(
	double a, double b,  /// Границы
	double h, /// Шаг
	int save_count, /// Количество элементов, которые надо сохрать в массив
	T ya, /// y(a)=ya
	std::function<T(double, const T&)> f, /// Функция f(t, y)
	bool is_write_progress = false /// Показывать ли прогресс расчета
);

//-----------------------------------------------------------------------------
/** Неявный метод Эйлера для решения обыкновенных ДУ: y'=f(t, y). */
template<class T>
std::vector<std::pair<double, T>> solveDE_Euler_Nonexplicit(
	double a, double b,  /// Границы
	double h, /// Шаг
	int save_count, /// Количество элементов, которые надо сохрать в массив
	T ya, /// y(a)=ya
	std::function<T(double, const T&)> f, /// Функция f(t, y)
	bool is_write_progress = false /// Показывать ли прогресс расчета
);

//-----------------------------------------------------------------------------
/** Метод Рунге-Кутта 4 порядка для решения обыкновенных ДУ: y'=f(t, y). */
template<class T>
std::vector<std::pair<double, T>> solveDE_Runge_Kutta4(
	double a, double b,  /// Границы
	double h, /// Шаг
	int save_count, /// Количество элементов, которые надо сохрать в массив
	T ya, /// y(a)=ya
	std::function<T(double, const T&)> f, /// Функция f(t, y)
	bool is_write_progress = false /// Показывать ли прогресс расчета
);

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
template<class T>
std::vector<std::pair<double, T>> solveDE_Euler_Explicit(double a, double b, double h, int save_count, T ya, std::function<T(double, const T&)> f, bool is_write_progress) {
    int n = (b-a)/h;

	std::vector<std::pair<double, T>> result;
	result.reserve(save_count);

    double t_n = a;
	T y_n = ya;
    for (int i = 0; i < n; ++i) {
		y_n = y_n + h*f(t_n, y_n);
		t_n += h;

		if (i % (n / save_count) == 0)
			result.push_back({t_n, y_n});

		if (is_write_progress && i % (n / 1000) == 0)
			std::cout << "\r" << std::setprecision(1) << std::fixed << std::setw(6) << i*100/double(n) << "%";
    }

    return result;
}

//-----------------------------------------------------------------------------
template<class T>
std::vector<std::pair<double, T>> solveDE_Euler_Nonexplicit(double a, double b, double h, int save_count, T ya, std::function<T(double, const T&)> f, bool is_write_progress) {
	int n = (b-a)/h;

	std::vector<std::pair<double, T>> result;
	result.reserve(save_count);

	double t_n = a;
	T y_n = ya;
	for (int i = 0; i < n; ++i) {
		double t_n1 = t_n + h;
		y_n = y_n + h/2.0 * (f(t_n, y_n) + f(t_n1, y_n + h*f(t_n, y_n)));
		t_n = t_n1;

		if (i % (n / save_count) == 0)
			result.push_back({t_n, y_n});

		if (is_write_progress && i % (n / 1000) == 0)
			std::cout << "\r" << std::setprecision(1) << std::fixed << std::setw(6) << i*100/double(n) << "%";
	}

	return result;
}

//-----------------------------------------------------------------------------
template<class T>
std::vector<std::pair<double, T>> solveDE_Runge_Kutta4(double a, double b, double h, int save_count, T ya, std::function<T(double, const T&)> f, bool is_write_progress) {
	int n = (b-a)/h;

	std::vector<std::pair<double, T>> result;
	result.reserve(save_count);

	double t_n = a;
	T y_n = ya;
	for (int i = 0; i < n; ++i) {
		double t_n1 = t_n + h;
		T k_n_1 = f(t_n, y_n);
		T k_n_2 = f(t_n + h / 2.0, y_n + h / 2.0 * k_n_1);
		T k_n_3 = f(t_n + h / 2.0, y_n + h /2.0 * k_n_2);
		T k_n_4 = f(t_n + h, y_n + h * k_n_3);

		y_n = y_n + h / 6.0 * (k_n_1 + 2.0 * k_n_2 + 2.0 * k_n_3 + k_n_4);
		t_n = t_n1;

		if (i % (n / save_count) == 0)
			result.push_back({t_n, y_n});

		if (is_write_progress && i % (n / 1000) == 0)
			std::cout << "\r" << std::setprecision(1) << std::fixed << std::setw(6) << double(i)*100/double(n) << "%";
	}

	return result;
}