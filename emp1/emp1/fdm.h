#pragma once
#ifndef FDM_H
#define FDM_H

#include "diagonal_matrix.h"
#include "mesh_generator.h"

class fdm
{
public:
	void mesh_to_slae(std::vector<node>& mesh);

	void calculate();

private:
	diagonal_matrix* A;		// Матрица СЛАУ
	std::vector<double> b;	// Вектор правой части
	std::vector<double> q;	// Результат

	uint32_t slae_size;		// Размер СЛАУ
};

#endif