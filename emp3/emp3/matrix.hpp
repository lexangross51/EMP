#pragma once
#ifndef MATRIX_H
#define MATRIX_H

#include "head.hpp"

enum class matrix_type
{
	DENSE,
	PROFILE,
	SPARSE
};

class matrix
{
	friend class solver;
	friend class mfe;

public:
	matrix(uint32_t _size) { dim = _size; }

	std::unique_ptr<std::vector<double>> dot(const std::vector<double>& vector);

	inline uint32_t size() const { return dim; };

	void to_dense();
	void to_sparse();
	void to_profile();

	void save(matrix_type type, std::string path = directory + "matrix\\");

private:
	uint32_t dim;

	std::vector<uint32_t> ig, jg;
	std::vector<double> di, ggl, ggu;

	std::vector<std::vector<double>> dense_matrix;
};

#endif