#define _USE_MATH_DEFINES

#include <fstream>
#include <iomanip>
#include <cmath>
#include <vector>

int main() {
	// Параметры метода
	double a = 0;
	double b = 100000;
	double h = 0.01; 
	int n = (b-a)/h;

	// Данные константы
	double qi = 1.0e-3;
	double pj = 1e5;

	// Константы, которые еще надо настроить
	double si = 1;
	double sj = si;
	double w = 1;
	double k = 1;
	double m = 1;

	// Функции производных
	auto vd = [si, pj, sj, w, m] (double pi, double v) -> double { return (pi*si - pj*sj - w*v)/m; };
	auto pid = [qi, si, k] (double v) -> double { return (qi - si*v)/k; };
	auto xd = [] (double v) -> double { return v; };

	// Массивы значений итоговых функций
	std::vector<double> v(n+1, 0);
	std::vector<double> pi(n+1, 0);
	std::vector<double> x(n+1, 0);

	// Начальные присвоения
	v[0] = 0;
	pi[0] = 1e5;
	x[0] = 0;

	// Итерации
	double t_n = a;
    for (int i = 1; i < n+1; ++i) {
        v[i] = v[i-1] + h*vd(pi[i-1], v[i-1]);
        pi[i] = pi[i-1] + h*pid(v[i-1]);
        x[i] = x[i-1] + h*xd(v[i-1]);

		t_n += h;
    }

	// Выводим результат
	std::ofstream fout("out.txt");
	fout << std::setprecision(16);
	for (size_t i = 0; i < v.size(); i+=1000)
		fout << v[i] << "\t" << pi[i] << "\t" << x[i] << std::endl;
	fout.close();
}