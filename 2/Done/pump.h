#pragma once

#include <vector>
#include <functional>
#include <fstream>
#include <iomanip>
#include <cmath>

using std::vector;


class pump {

public:
	pump();
	double speed(double Pi, double V) { return (Pi*Si - Pj*Sj - nu*V) / m; }
	double xCoordinate(double V) { return V; }
	double pressure(double V) { return (qi - Si*V) / K_elast; }

	// Выводим вектора X, V, Pi
	void outputResults(char *fileName);

	// Явный метод Эйлера для решения обыкновенных ДУ
	void solveDE_Euler_Explicit();

	// Неявный метод Эйлера для решения обыкновенных ДУ
	void solveDE_Euler_Nonexplicit();

	// Метод Рунге-Кутта 4 порядка для решения обыкновенных ДУ
	void solveDE_Runge_Kutta4();


private:
	// Параметры метода
	double a = 0;				// левый ограничитель
	double b = 10;				// правый ограничитель
	double h = 0.01;			// шаг
	int n = (b - a) / h;		// число шагов

	// Данные константы
	const double qi = 1.0e-3;	// расход воды слева
	const double Pj = 1e5;		// (атмосферное) давление справа

	// Константы, которые еще надо настроить
	const double Si = 1;		// площадь цилиндра слева
	const double Sj = 1;		// площадь цилиндра справа
	const double w = 1;			// 
	const double K_elast = 1;	// коэффициент упругости
	const double m = 1;			// масса поршня
	const double nu = 1;		// ню

	// Массивы значений итоговых функций
	vector <double> X, V, Pi;
};