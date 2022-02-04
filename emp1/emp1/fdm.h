#pragma once
#ifndef FDM_H
#define FDM_H

#include "diagonal_matrix.h"
#include "mesh.h"
#include "solver.h"
#include "structures.h"

class fdm
{
public:
	void mesh_to_slae(mesh& mesh, func2D u, function2D f);

	void calculate();

private:
	diagonal_matrix* A;		// Матрица СЛАУ
	std::vector<double> b;	// Вектор правой части
	std::vector<double> q;	// Результат

	double beta = 0;
	double u_beta = 0;

	uint32_t slae_size;		// Размер СЛАУ

	void write(	const std::string dir, 
				const double w, 
				const std::pair<uint32_t, double> result);

	double du_dx(func2D& f, point& p);
	double du_dy(func2D& f, point& p);
};

#endif