#pragma once
#ifndef MFE_HPP
#define MFE_HPP

#include "mesh.hpp"
#include "matrix.hpp"
#include "basis_function.hpp"
#include "gauss.hpp"
#include "portrait_generator.hpp"

class mfe
{
public:
	mfe(space_grid& mesh);

	void assembly_global_matrix_and_vector(r_function3D& fs, r_function3D& fc);

	inline void set_w(const double _w) { w = _w; }

private:
	double w;

	std::vector<std::vector<double>> local_p;
	std::vector<std::vector<double>> local_c;

	std::vector<double> local_vec_fs;
	std::vector<double> local_vec_fc;

	void build_local_matrices(finite_elem& elem);
	void build_local_vector(finite_elem& elem, r_function3D& fs, r_function3D& fc);

	void add_to_global_matrix(const uint32_t i, const uint32_t j, const double val);

	space_grid* grid;		// Сетка

	matrix* A;				// Глобальная матрица
	std::vector<double> q;	// Вектор весов
	std::vector<double> b;	// Вектор правой части

	gauss gauss;			// Гаусс (для интегрирования)
};

#endif