#pragma once
#include <vector>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <iomanip>

const std::string directory = "C:\\Users\\lexan\\OneDrive\\Рабочий стол\\НГТУ\\6 семестр\\Уравнения математической физики\\emp1\\emp1\\";

struct point
{
	double x;
	double y;

	point(double _x = 0.0, double _y = 0.0) :
		x(_x), y(_y) {};
};

// Граница: задает пределы (2 точки) и тип краевого
struct border
{
	enum class bound_cond
	{
		NONE,
		DIRICHLET,
		NEUMANN,
		NEWTON
	};

	point limits[2];
	bound_cond bc = bound_cond::NONE;
};

// Элементарная область: 4 границы и коэффициенты
struct area
{
	border borders[4];
	double lambda = 0.0;
	double gamma = 0.0;
};

// Узел сетки: 
// координаты, 
// позиции в массивах X и Y, 
// тип (внутренний, граничный, фиктивный)
// тип краевого
// коэффициенты лямбда и бета
struct node
{
	enum class node_type 
	{ 
		INTERNAL, 
		BORDER, 
		FICTITIOUS 
	};

	point p;
	uint32_t x_pos, y_pos;
	node_type type;
	border::bound_cond bc;
	double lambda;
	double gamma;

	node(	
			point& _p, uint32_t i, uint32_t j, 
			node_type _type = node_type::FICTITIOUS, 
			border::bound_cond _bc = border::bound_cond::NONE,
			double _lambda = 0.0, double _gamma = 0.0
		) :
		p(_p), x_pos(i), y_pos(j), type(_type), bc(_bc), lambda(_lambda), gamma(_gamma) {};
};