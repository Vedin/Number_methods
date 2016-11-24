#ifndef FILE_OPS
#define FILE_OPS

#include <vector>
#include <fstream>
using namespace std;

#define Vector vector<double>

void read_dots_from_file(const char* filename, Vector &x, Vector &fx) {
	x.clear();
	fx.clear();

	ifstream in(filename);

	// Читаємо кількість вузлів
	int n; in >> n;
	
	double tmp;

	// Читаємо Х
	for (int i = 0; i < n; i++) {
		in >> tmp; x.push_back(tmp);
	}

	// Читаємо f(X)
	for (int i = 0; i < n; i++) {
		in >> tmp; fx.push_back(tmp);
	}
}
void create_chebyshev_dots_file( double (*f)(double), int n, const char* filename) {
	const double pi = 3.1415926;

	ofstream out(filename);
	Vector fx;

	out << n << endl;

	for (int i = 1; i <= n; i++) {
		
		// Рахуєм корені многочлена Чебишева
		double x = cos(((2 * i) - 1) * pi / (2 * n));
		
		// Зберігаєм результат функції від цієї точки
		fx.push_back(f(x));
		
		// Виводимо
		out << x << " ";
	}
	out << endl;

	for (int i = 0; i < n; i++) {
		out << fx.at(i) << " ";
	}

}
void create_simple_dots_file( double (*f)(double), int n, const char* filename, double hi, double lo) {
	const double pi = 3.1415926;

	double step = (hi - lo) / (n + 1);

	ofstream out(filename);
	Vector fx;

	out << n << endl;

	int toGen = n;

	for (int i = 1; i < (n + 1); i++) {
		double x = lo + i * step;

		if (x == 0) {
			++n;
			continue;
		}
		fx.push_back(f(x));

		out << x << " ";
	}

	out << endl;

	for (int i = 0; i < toGen; i++) {
		out << fx.at(i) << " ";
	}

}

#undef Vector

#endif