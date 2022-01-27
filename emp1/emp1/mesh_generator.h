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
	uint32_t n_omega;			// Количество подобластей
	uint32_t nx, ny;			// Кол-во элементов в массивах 
								// X_lines и Y_lines соответственно
	uint32_t num_x, num_y;		// Число по разбиений по X и Y 
								// соответственно
	double kx, ky;				// Коэф разрядки по X и Y 
								// соответственно
	std::vector<area> areas;	// Массив с подобластями
		
	std::vector<double> X_lines;	// Координатные линии по X
	std::vector<double> Y_lines;	// Координатные линии по Y
	std::vector<double> X;			// Массив X с учетом разбиений
	std::vector<double> Y;			// Массив Y с учетом разбиений

	void input(const std::string dir);

	void generate_xy();

	node::node_type what_type(point &p);

	bool is_in_area(const point& p, uint32_t &area_num);

	void write_by_type(const std::vector<node>& mesh, const std::string dir);
};

#endif