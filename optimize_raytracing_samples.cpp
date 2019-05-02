/**
	Программа для нахождения максимально равномерно расположенных точек внутри квадрата. Будет использоваться для рейтрейсинга с задаваемым числом семплов. Число семплов от 1 до 512.

	Для этого использовалась идея что точки можно снабдить отталкивающей гравитацией и запустить симуляцию, что это может сойтись к стабильной конфигурации.

	Результат: точки осциллируют и ни к чему не сходятся. Результат показан на гифке gravity.gif.

 */

#include <twg/image/image_drawing.h>

#undef min
#undef max

#include <vector>
#include <random>
#include <fstream>
#include <cmath>
#include "3/methods.h"
#include "3/vec.h"

#include <spob/spob.h>

using namespace std;

//-----------------------------------------------------------------------------
double random(void) {
  static mt19937 generator;
  static uniform_real_distribution<double> distribution(0, 1);
  return distribution(generator);
}

//-----------------------------------------------------------------------------
vec f(double t, const vec& in) {
	const double restriction_force = 100;
	const double stability_eps = 0.001;

	int n = in.size() / 4;

	vector<vec> v, x;
	for (int i = 0; i < 4*n; i += 4) {
		v.push_back({in[i+0], in[i+1]});
		x.push_back({in[i+2], in[i+3]});
	}

	vector<int> x_neigh = {-1, 0, 1, 1, 1, 0, -1, -1};
	vector<int> y_neigh = {-1, -1, -1, 0, 1, 1, 1, 0};

	vector<vec> neighbors;
	for (int i = 0; i < 8; ++i) {
		for (auto& j : x) {
			neighbors.push_back({j[0] + x_neigh[i], j[1] + y_neigh[i]});
		}
	}

	x.insert(x.end(), neighbors.begin(), neighbors.end());

	vec result(in.size(), 0);
	for (int i = 0; i < 4*n; i += 4) {
		auto px = in[i+2];
		auto py = in[i+3];
		auto& vx = result[i+0];
		auto& vy = result[i+1];
		for (auto& j : x) {
			#define sqr(a) ((a)*(a))
			double len = sqr(px-j[0]) + sqr(py-j[1]);
			#undef sqr
			if (len != 0) {
				vx += (j[0]-px)/len;
				vy += (j[1]-py)/len;
			}
		}

		// Делаем из притяжения отталкивание
		vx = -vx;
		vy = -vy;

		/*if (std::fabs(vx + vy) < stability_eps) {
			vx = ((in[i+0] > 0) ? -1 : 1) * restriction_force;
			vy = ((in[i+1] > 0) ? -1 : 1) * restriction_force;
		}*/

		// Чтобы не выходило за границы
		if (px > 1) 
			vx = -restriction_force; 
		if (px < 0) 
			vx = restriction_force;
		if (py > 1) 
			vy = -restriction_force;
		if (py < 0) 
			vy = restriction_force;

		result[i+2] = in[i+0];
		result[i+3] = in[i+1];
	}

	return result;
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

void drawCircle(twg::ImageDrawing_aa& img, twg::Point_d pos, double r) {
	img.drawPolygon(twg::computeEllipse({r, r}).move(pos));
}

int main() {
	int size = 50;

	vec start(4*size);
	for (int j = 0; j < 4*size; j += 4) {
		start[j+0] = 0;
		start[j+1] = 0;
		start[j+2] = 0.5 +
			0.25 * sin(2.0 * 3.14159265359 / double(size) * double(j/4)) +
			0.1 * random();
		start[j+3] = 0.5 +
			0.25 * cos(2.0 * 3.14159265359 / double(size) * double(j/4)) + 
			0.1 * random();
	}
	double a = 0;
	double b = 1;
	double n = 5000;
	double h = (b - a) / n;
	auto result = solveDE_Runge_Kutta4<vec>(a, b, h, 100, start, f, true);

	using namespace twg;

	ImageGif gif;

	gif.start({500, 500}, "gravity.gif");

	ImageDrawing_aa img({500, 500});
	img.setBrush(setAlpha(Black, 128));

	img.clear(White);

	spob::space2 space({3, 0}, {0, 3}, {-1, -1});
	spob::space2 screen({500, 0}, {0, 500}, {0, 0});

	spob::space2 toScreen = combine(invert(screen), space);

	int counter = 0;
	for (auto& i : result) {
		img.clear();
		img.drawText({0, 0}, std::to_wstring(counter++));

		auto draw_line = [&] (spob::vec2 a, spob::vec2 b) {
			img.drawLine(screen.from(space.to(a)), screen.from(space.to(b)));
		};
		draw_line({-1, 0}, {2, 0});
		draw_line({-1, 1}, {2, 1});

		draw_line({0, -1}, {0, 2});
		draw_line({1, -1}, {1, 2});

		auto in = i.second;

		vector<vec> v, x;
		for (int i = 0; i < 4*size; i += 4) {
			v.push_back({in[i+0], in[i+1]});
			x.push_back({in[i+2], in[i+3]});
		}

		vector<int> x_neigh = {-1, 0, 1, 1, 1, 0, -1, -1};
		vector<int> y_neigh = {-1, -1, -1, 0, 1, 1, 1, 0};

		vector<vec> neighbors;
		for (int i = 0; i < 8; ++i) {
			for (auto& j : x) {
				neighbors.push_back({j[0] + x_neigh[i], j[1] + y_neigh[i]});
			}
		}

		x.insert(x.end(), neighbors.begin(), neighbors.end());

		for (auto& j : x) {
			spob::vec2 pos(j[0], j[1]);
			pos = space.to(pos);
			pos = screen.from(pos);
			drawCircle(img, pos, 3);
		}

		/*ofstream fout("samples.txt");
		fout << size << endl;
		for (auto& i : x) {
			fout << "(" << i[0] << ", " << i[1] << ")" << endl;
		}
		fout << endl;
		fout.close();*/

		gif.process(img, 2);
	}

	gif.end();
}