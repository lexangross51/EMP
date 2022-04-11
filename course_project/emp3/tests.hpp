#pragma once
#ifndef TESTS_HPP
#define TESTS_HPP

#include "mesh.hpp"
#include "mesh_generator.hpp"
#include "mfe.hpp"

class tests
{
public:
    tests();

    void run();

private:
    std::vector<function2D> _u_rz;
    std::vector<function1D> _u_t;
    std::vector<std::vector<function2D_t_r>> _f;

    std::vector<std::string> _u_rz_names;
    std::vector<std::string> _u_t_names;
};

#endif