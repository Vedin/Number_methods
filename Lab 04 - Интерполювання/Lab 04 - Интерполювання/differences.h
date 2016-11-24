#ifndef DIFF_H
#define DIFF_H

#include <vector>
#include <iostream>
#include <fstream>
using namespace std;

#define Matrix vector<vector<double>>
#define Vector vector<double>

// Чи всі числа у векторі однакові
bool all_numbers_are_same(Vector nums, int lastpos) {

	for (int i = 0; i < lastpos - 1; i++) {
		if (nums.at(i) != nums.at(i + 1)) 
			return false;
	}

	return true;
}
Matrix build_diff_matrix (Vector x, Vector fx) {
	Matrix result;			// Зберігає результат
	result.push_back(fx);	// В перший рядок перепишемо наші значення функції
	int max_size = fx.size(); 


	int line_to_form = 0;
	Vector temp_nums;

	// Формуємо наступні рядки
	do {
		++line_to_form;
		temp_nums.clear();

		for (int i = 0; i < max_size - line_to_form; i++) {
			int first_x_index = i; 
			int second_x_index = i + line_to_form;

			double diff = result.at(line_to_form - 1).at(i + 1) -
						  result.at(line_to_form - 1).at(i);

			diff /= (x.at(second_x_index) - x.at(first_x_index));

			temp_nums.push_back(diff);
		}

		result.push_back(temp_nums);
	} 
	// Поки не знайдемо рядок, в якому всі числа однакові
	
	while (line_to_form != max_size - 1);
	
	//while (all_numbers_are_same(temp_nums, max_size - line_to_form) == false);
	
	return result;
}

double get_diff(const Matrix & diffs, int first_i, int second_i) {
	if (first_i > second_i) swap(first_i, second_i);

	return diffs.at(second_i - first_i).at(first_i);
}


#undef Matrix
#undef Vector

#endif