#include "head.hpp"

struct finite_element
{
	std::vector<uint32_t> nodes;

	double lambda, sigma, hi;

	uint32_t operator[](const uint32_t index) { return nodes[index]; }

	finite_element(std::vector<uint32_t> _nodes = {}, double _lambda = 0.0, double _sigma = 0.0, double _hi = 0.0)
	{
		nodes = _nodes;
		lambda = _lambda;
		sigma = _sigma;
		hi = _hi;
	}
};