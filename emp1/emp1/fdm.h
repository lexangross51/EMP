#pragma once
#ifndef FDM_H
#define FDM_H

#include "diagonal_matrix.h"
#include "mesh.h"
#include "solver.h"

class fdm
{
public:
	void mesh_to_slae(mesh& mesh, func2D_u& u, func2D_f& f);

	std::pair<uint32_t, double> calculate(std::vector<double>& q);

private:
	diagonal_matrix* A;			// Матрица СЛАУ
	std::vector<double> b;		// Вектор правой части
	std::vector<double> q;		// Результат

	std::vector<double> exact;	// Точное решение

	double beta = 0;
	double u_beta = 0;

	uint32_t slae_size;			// Размер СЛАУ

	double du_dx(func2D_u& f, point& p);
	double du_dy(func2D_u& f, point& p);
};

#endif