#pragma once
#ifndef SOLVER_H
#define SOLVER_H

#include "matrix.h"

enum iterative_method
{
	LOS_LU,
	BCGSTAB_LU
};

class solver
{
public:
	std::pair<uint32_t, double> solve_iterative(
		matrix& A, 
		std::vector<double>& b, 
		std::vector<double>& x, 
		iterative_method method,
		uint32_t max_iter, double eps
	);

	void solve_by_LU(matrix& A, std::vector<double>& b, std::vector<double>& x);

private:
	iterative_method method;

	// Параметры для итерационного решателя
	uint32_t max_iter;
	double eps;

	// Итерационные методы
	std::pair<uint32_t, double> LOS_LU(matrix& A,
		std::vector<double>& b,
		std::vector<double>& x,
		uint32_t max_iter, double eps
	);

	std::pair<uint32_t, double> BCGSTAB_LU(matrix& A,
		std::vector<double>& b,
		std::vector<double>& x,
		uint32_t max_iter, double eps
	);
};

#endif