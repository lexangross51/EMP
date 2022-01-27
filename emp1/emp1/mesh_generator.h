#pragma once
#ifndef MESH_GENERATOR_H
#define MESH_GENERATOR_H

#include "head.h"

class mesh_generator
{
public:
	mesh_generator()
	{
		n_omega = 0;
		nx = ny = 0;
		num_x = num_y = 0;
		kx = ky = 0.0;
	}

	void build_mesh(std::vector<node>& mesh);

private:
	uint32_t n_omega;			// ���������� �����������
	uint32_t nx, ny;			// ���-�� ��������� � �������� 
								// X_lines � Y_lines ��������������
	uint32_t num_x, num_y;		// ����� �� ��������� �� X � Y 
								// ��������������
	double kx, ky;				// ���� �������� �� X � Y 
								// ��������������
	std::vector<area> areas;	// ������ � ������������
		
	std::vector<double> X_lines;	// ������������ ����� �� X
	std::vector<double> Y_lines;	// ������������ ����� �� Y
	std::vector<double> X;			// ������ X � ������ ���������
	std::vector<double> Y;			// ������ Y � ������ ���������

	void input(const std::string dir);

	void generate_xy();

	node::node_type what_type(point &p);

	bool is_in_area(const point& p, uint32_t &area_num);

	void write_by_type(const std::vector<node>& mesh, const std::string dir);
};

#endif