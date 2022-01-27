#pragma once
#include <vector>
#include <fstream>
#include <iostream>
#include <algorithm>

const std::string directory = "C:\\Users\\lexan\\source\\repos\\emp1\\emp1\\";

struct point
{
	double x;
	double y;

	point(double _x = 0.0, double _y = 0.0) :
		x(_x), y(_y) {};
};

//  --2
// |  |
// |  |
// 1--
struct area
{
	point borders[2];
	double lambda = 0.0;
	double gamma = 0.0;
};

struct node
{
	enum class node_type 
	{ 
		INTERNAL, 
		BORDER, 
		FICTITIOUS 
	};

	point p;
	node_type type;
	double lambda;
	double gamma;

	node(point& _p, node_type _type = node_type::FICTITIOUS, 
		 double _lambda = 0.0, double _gamma = 0.0) :
		p(_p), type(_type), lambda(_lambda), gamma(_gamma) {};
};
