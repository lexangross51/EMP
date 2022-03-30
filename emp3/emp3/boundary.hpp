#pragma once
#ifndef BOUNDARY_HPP
#define BOUNDARY_HPP

#include "head.hpp"

enum class boundary_type {
	DIRICHLET = 1,
	NEUMANN = 2
};

struct boundary_cond {
	uint8_t edge;
	boundary_type type;
};

struct dirichlet_cond {
	uint32_t node;
	double value;
};

struct neumann_cond {
	uint32_t ielem;
	uint32_t local_node_1, local_node_2;
	double tetta;
};

#endif