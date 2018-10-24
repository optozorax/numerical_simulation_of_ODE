#define _USE_MATH_DEFINES

#include <fstream>
#include <iomanip>
#include <cmath>
#include <vector>

int main() {
	// Параметры метода
	double a = 0;
	double b = 100;
	double h = 0.0001; 
	int n = (b-a)/h;

	// Данные константы
	const double qi = 1.0e-3;
	const double pj = 1e5;

	// Константы, которые еще надо настроить
	const double si = 1;
	const double sj = 1;
	const double w = 1;
	const double k = 1;
	const double m = 1;

	// Массивы значений итоговых функций
	std::vector<double> v(n+1, 0);
	std::vector<double> pi(n+1, 0);
	std::vector<double> x(n+1, 0);

	// Начальные значения
	v[0] = 0;
	pi[0] = 1e5;
	x[0] = 0;

	// Итерации
	double t_n = a;
    for (int i = 1; i < n+1; ++i) {
		const double t0 = t_n;
		const double v0 = v[i-1];
		const double x0 = x[i-1];
		const double t1 = t_n + h;
		double& v1 = v[i];
		double& x1 = x[i];
		double& p1 = pi[i];

		const double a = qi * si / (m * k);
		const double b = si * si / (m * k);
		const double c = w / m;
		const double d = (pj * sj - 1e5 * si) / m;
		const double y = t0*a - x0*b - v0*c - d;

		auto euler_1 = [&] () {
			v1 = v0 + h*y;
			x1 = x0 + h*v0;
		};

		auto euler_2 = [&] () {
			v1 = -(2.0*y*h*(c*h - 1.0) + 2.0*h*(b*x0 + d) + v0*(b*h*h + 2.0*c*h - 4.0) - 2.0*a*h*t1)/(b*h*h + 4.0);
			x1 = x0 + h/2.0*(v0+v1);
		};

		euler_1();
		//euler_2();

		p1 = (qi*t1 - si*x1)/k + pi[0];

		t_n += h;
    }

    // Удаляем начальные значения
	// v.erase(v.begin());
	// pi.erase(pi.begin());
	// x.erase(x.begin());

	// Выводим результат
	std::ofstream fout("out.txt");
	fout << std::setprecision(16);
	for (size_t i = 0; i < v.size(); i+=1000)
		fout << v[i] << "\t" << pi[i] << "\t" << x[i] << std::endl;
	fout.close();
}