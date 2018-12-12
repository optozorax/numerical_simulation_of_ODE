#pragma once

// Данные константы
const double q_nasos = 1.0e-3;	// расход воды слева
const double P_atm = 1e5;		// (атмосферное) давление справа
const double pi = 3.141592653589793238463;

// Константы, которые еще надо настроить
const double s_i = 1; 		// площадь цилиндра слева
const double s_j = 1;	// площадь цилиндра справа
const double w = 1;			// 
const double K_elast = 1e-5;	// коэффициент упругости
const double m = 1;			// масса поршня
const double nu = 0.01;		// ню
const double x_swap_r = 0.009;
const double x_stop_r = x_swap_r + 0.001;
const double x_swap_l = -0.001;
const double x_stop_l = x_swap_l - 0.001;

const double d = 2e-3;		//
const double l = 5e-3;		//
const double v = 10;// l*s_i;		//
const double gamma = 1.4;	//
const double pAtm = 1e5;	// атмосферное давление
const double ro = 1e3;		// плотность воды
const double ReCnt = 321;
const double Nu = 1e-6;		// число Нуссельта воды
const double dH = pi * d / (4 * s_i);
const double kappa = 1;

const double Es = 1260 * 1260 * 1000; 		// модуль жёсткость воды ~10^9
const double C = v / Es;					// жёсткость воды
const double C_cav = pow(pAtm, 1 / gamma);	// жёсткость воды при кавитации
const double mltRe = dH / Nu * s_i;
const double r = 12 * ro*Nu*l / (dH*dH*s_i);
const double r_kappa = 0.1582 * ro * pow(Nu, 0.25)*l / (pow(dH, 1.25)*pow(s_i, 1.75));

const double B = 1 / (l*sqrt(2 * ro));
const double F = s_i * sqrt(2 * ro);
const double eps = 0.57 + 0.043*(1.1 - 1);
const double dzeta = 0.035*(pow(1 / eps - 1, 2));