#pragma once
#ifndef FDM_H
#define FDM_H

#include "diagonal_matrix.h"
#include "mesh.h"
#include "solver.h"

class fdm
{
public:
	void mesh_to_slae(mesh& mesh, function2D u, function2D f);

	void calculate();

private:
	diagonal_matrix* A;		// Матрица СЛАУ
	std::vector<double> b;	// Вектор правой части
	std::vector<double> q;	// Результат

	uint32_t slae_size;		// Размер СЛАУ

	void write(	const std::string dir, 
				const double w, 
				const std::pair<uint32_t, double> result);
};

#endif