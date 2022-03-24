#pragma once
#ifndef MESH_H
#define MESH_H

#include "head.h"

struct point2D
{
	double x, y;

	point2D(double _x, double _y) :
		x(_x), y(_y) {};
};

struct finite_elem
{
	uint32_t number;
	std::vector<uint32_t> nodes;

	double lambda, sigma, hi;

	uint32_t operator[](const uint32_t index) { return nodes[index]; }

	finite_elem(uint32_t _number, std::vector<uint32_t> _nodes, double _lambda, double _sigma, double _hi)
	{
		number = _number;
		nodes = _nodes;
		lambda = _lambda;
		sigma = _sigma;
		hi = _hi;
	}
};

class mesh
{
public:
	inline const finite_elem& operator[](const uint32_t index) { return elems[index]; }

	inline point2D& point(const uint32_t index) { return points[index]; }

	inline uint32_t size(void) const { return dimension; }

	void save(std::string path = directory + "mesh\\");

private:
	uint32_t dimension;
	std::vector<finite_elem> elems;
	std::vector<point2D> points;
};

#endif