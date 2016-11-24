#include "graph_plotter.h"
#include "file_operations.h"
#include "polinom.h"

#define Vector vector<double>
#define Matrix vector<Vector>

double f(double x) {
	return log(cos(x))/log(2.);
}

const double pi = 3.1415926;


int main() {

	Vector x;
	Vector fx;

	double global_scale = 50;

	create_chebyshev_dots_file(f, 3, "chebyshev_dots.txt");
	create_simple_dots_file(f, 3, "simple_dots.txt", pi / 2, - pi / 2);

	read_dots_from_file("chebyshev_dots.txt", x, fx);
	Polinom chebyshev(x, fx);

	read_dots_from_file("simple_dots.txt", x, fx);
	Polinom simple(x, fx);

	auto cheb_res = plotGraph(chebyshev, pi/2, - pi/2, global_scale);
	auto simp_res = plotGraph(simple, pi/2, - pi/2, global_scale);
	auto orig_res = plotGraph(f, pi/2, - pi/2, global_scale);

	auto console = GetConsoleWindow();
	HDC handle = GetDC(console);

	int zero_y_level = 100;
	int zero_x_level = 100;

	drawGraph(handle, cheb_res, global_scale, GREEN, zero_x_level, zero_y_level);
	drawGraph(handle, simp_res, global_scale, RED, zero_x_level + 200, zero_y_level);
	drawGraph(handle, orig_res, global_scale, BLUE, zero_x_level + 400, zero_y_level);

	system("Pause");
	return 0;
}