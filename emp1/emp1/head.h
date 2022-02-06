#pragma once
#include <vector>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <functional>
#include <string>

const std::string directory = "C:\\Users\\lexan\\OneDrive\\Рабочий стол\\НГТУ\\6 семестр\\Уравнения математической физики\\emp1\\emp1\\";

typedef std::function<double(double, double)> func2D_u;
typedef std::function<double(double, double, double, double)> func2D_f;

struct point
{
	double x;
	double y;

	point(double _x = 0.0, double _y = 0.0) :
		x(_x), y(_y) {};
};