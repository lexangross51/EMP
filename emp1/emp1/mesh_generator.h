#pragma once
#ifndef MESH_GENERATOR_H
#define MESH_GENERATOR_H

#include "structures.h"

class mesh_generator
{
public:
	enum class grid_type
	{
		UNIFORM,
		NONUNIFORM
	};

	mesh_generator()
	{
		n_omega = 0;
		nx = ny = 0;
	}

	void build_mesh(std::vector<node>& mesh, grid_type type);

private:
	uint32_t n_omega;						 // ���-�� �����������
	uint32_t nx, ny;						 // ���-�� ��������� � �������� 
											 // X_lines � Y_lines ��������������

	std::vector<area> areas;				 // ������ � ������������
		
	std::vector<double> X_lines;			 // ������������ ����� �� X
	std::vector<double> Y_lines;			 // ������������ ����� �� Y

	std::vector<std::vector<uint32_t>> part; // ������ � ����������� � ����������
	std::vector<std::vector<double>> kr;	 // ������ � �������������� ��������

	std::vector<double> X;					 // ������ X � ������ ���������
	std::vector<double> Y;					 // ������ Y � ������ ���������

	void input(const std::string dir);

	void generate_xy(const grid_type type);

	node::node_type what_type(const point &p, const uint32_t i, const uint32_t j);

	uint32_t what_border(const point& p, const uint32_t area_num);

	bool is_in_area(const point& p, uint32_t &area_num);

	void write_by_type(const std::vector<node>& mesh, const std::string dir);
};

#endif