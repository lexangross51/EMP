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
	uint32_t n_omega;						 // Кол-во подобластей
	uint32_t nx, ny;						 // Кол-во элементов в массивах 
											 // X_lines и Y_lines соответственно

	std::vector<area> areas;				 // Массив с подобластями
		
	std::vector<double> X_lines;			 // Координатные линии по X
	std::vector<double> Y_lines;			 // Координатные линии по Y

	std::vector<std::vector<uint32_t>> part; // Массив с информацией о разбиениях
	std::vector<std::vector<double>> kr;	 // Массив с коэффициентами разрядки

	std::vector<double> X;					 // Массив X с учетом разбиений
	std::vector<double> Y;					 // Массив Y с учетом разбиений

	void input(const std::string dir);

	void generate_xy(const grid_type type);

	node::node_type what_type(const point &p, const uint32_t i, const uint32_t j);

	uint32_t what_border(const point& p, const uint32_t area_num);

	bool is_in_area(const point& p, uint32_t &area_num);

	void write_by_type(const std::vector<node>& mesh, const std::string dir);
};

#endif