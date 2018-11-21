#pragma once

#include <vector>
#include <functional>

//-----------------------------------------------------------------------------
/** Явный метод Эйлера для решения обыкновенных ДУ: y'=f(t, y). */
typedef<template T>
std::vector<T> solveDE_Euler_Explicit(
	double a, double b,  /// Границы
	double h, /// Шаг
	T ya, /// y(a)=ya
	std::function<T(double, T)> f /// Функция f(t, y)
);

//-----------------------------------------------------------------------------
/** Неявный метод Эйлера для решения обыкновенных ДУ: y'=f(t, y). */
typedef<template T>
std::vector<T> solveDE_Euler_Nonexplicit(
	double a, double b,  /// Границы
	double h, /// Шаг
	T ya, /// y(a)=ya
	std::function<T(double, T)> f /// Функция f(t, y)
);

//-----------------------------------------------------------------------------
/** Метод Рунге-Кутта 4 порядка для решения обыкновенных ДУ: y'=f(t, y). */
typedef<template T>
std::vector<T> solveDE_Runge_Kutta4(
	double a, double b,  /// Границы
	double h, /// Шаг
	T ya, /// y(a)=ya
	std::function<T(double, T)> f /// Функция f(t, y)
);