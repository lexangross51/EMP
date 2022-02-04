#pragma once
#include <vector>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <functional>

const std::string directory = "C:\\Users\\lexan\\OneDrive\\Рабочий стол\\НГТУ\\6 семестр\\Уравнения математической физики\\emp1\\emp1\\";

typedef std::function<double(double, double, double)> function2D;
typedef std::function<double(double, double)> func2D;

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