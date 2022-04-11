#pragma once
#ifndef BOUNDARY_HPP
#define BOUNDARY_HPP

#include "head.hpp"

enum class boundary_type {
    NONE,
    DIRICHLET,
    NEUMANN
};

struct boundary_cond {
    uint8_t border = 1;
    boundary_type type = boundary_type::NONE;
    dvector boundaries;
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