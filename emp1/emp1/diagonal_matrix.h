#pragma once
#ifndef DIAGONAL_MATRIX_H
#define DIAGONAL_MATRIX_H

#include "structures.h"

class diagonal_matrix
{
public:
	diagonal_matrix(const uint32_t _dim, const uint32_t _diags_count);

	inline uint32_t get_size(void) const { return dim; }

	inline double get_diag_elem(const uint32_t di, const uint32_t elem) { return diags[di][elem]; }

	std::unique_ptr<std::vector<double>> dot(const std::vector<double>& vector);

	double dot(const uint32_t row, const std::vector<double>& vector);

	void to_dense(std::string file_to_save);

public:
	uint32_t dim;								// Размер матрицы
	std::vector<std::vector<double>> diags;		// Диагонали
	std::vector<int> offset;					// Смещения побочных диагоналей
												// относительно главной
	std::vector<std::vector<double>> matrix;	// Матрица в плотном формате

	void make_dom();
};

#endif