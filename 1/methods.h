#pragma once

#include <vector>
#include <functional>
#include <iomanip>
const int _PRECISION_ = 5;


typedef std::function<double(double, double)> DE_F;
typedef std::function<double(double)> DE_Solve;

//-----------------------------------------------------------------------------
/** Явный метод Эйлера для решения обыкновенных ДУ: y'=f(t, y). */
std::vector<double> solveDE_Euler_Explicit(
	double a, double b,  /// Границы
	double h, /// Шаг
	double ya, /// y(a)=ya
	DE_F f /// Функция f(t, y)
);

//-----------------------------------------------------------------------------
/** Неявный метод Эйлера для решения обыкновенных ДУ: y'=f(t, y). */
std::vector<double> solveDE_Euler_Nonexplicit(
	double a, double b,  /// Границы
	double h, /// Шаг
	double ya, /// y(a)=ya
	DE_F f /// Функция f(t, y)
);

//-----------------------------------------------------------------------------
/** Метод Рунге-Кутта 4 порядка для решения обыкновенных ДУ: y'=f(t, y). */
std::vector<double> solveDE_Runge_Kutta4(
	double a, double b,  /// Границы
	double h, /// Шаг
	double ya, /// y(a)=ya
	DE_F f /// Функция f(t, y)
);