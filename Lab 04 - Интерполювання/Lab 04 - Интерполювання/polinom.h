#ifndef POLINOM_H
#define POLINOM_H

#include "differences.h"

#define Vector vector<double>
#define Matrix vector<Vector>

class Polinom {
	Vector x_values;
	Vector fx_values;

	Matrix diff;


	double get_polinom_val(int cur_x_index, double x) {

		if (cur_x_index >= x_values.size()) return 0;
		
		return get_diff(diff, 0, cur_x_index) +		// Розділена різниця +
			   (x - x_values.at(cur_x_index)) *		// (Х - Хі) *
			   get_polinom_val(cur_x_index + 1, x); // наступна дужка
	}

public:

	Polinom(Vector nx, Vector nfx) {
		x_values = nx;
		fx_values = nfx;

		diff = build_diff_matrix(x_values, fx_values);
	
	}

	double operator() (double x) {

		return get_polinom_val(0, x);
	}
};

#undef Matrix
#undef Vector

#endif