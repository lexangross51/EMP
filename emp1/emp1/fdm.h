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
	diagonal_matrix* A;		// ������� ����
	std::vector<double> b;	// ������ ������ �����
	std::vector<double> q;	// ���������

	uint32_t slae_size;		// ������ ����
};

#endif