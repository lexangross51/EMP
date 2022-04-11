#pragma once
#ifndef MFE_HPP
#define MFE_HPP

#include "basis_function.hpp"
#include "mesh.hpp"
#include "matrix.hpp"
#include "solver.hpp"
#include "gauss.hpp"
#include "portrait_generator.hpp"

class mfe
{
public:
    mfe(space_grid* space, time_grid* time, function2D_t& u, function2D_t_r& f);

    double solve();

private:
    space_grid* _mesh_rz;
    time_grid* _time_layers;

    uint8_t _local_size;
    std::vector<std::vector<double>> _local_G;
    std::vector<std::vector<double>> _local_M;
    std::vector<double> _local_b;

    function2D_t _u;
    function2D_t_r _f;

    matrix* _global_A;

    dvector _global_b;
    dvector _q;

    dvector _q_prev3;
    dvector _q_prev2;
    dvector _q_prev1;

    gauss _gauss;

    dvector _exact;

    void rebuild_mesh(space_grid* grid);

    void build_local_matrices(finite_element& elem);

    void build_local_vector(finite_element& elem, double time);

    void add_to_global_matrix(const uint32_t i, const uint32_t j, const double val);

    void initial_condition();

    void assembly_global_matrix_and_vector(const uint32_t time_m);

    void add_dirichlet(const uint32_t time_m);

    void add_meumann(const uint32_t time_m);

    void calc_exact(const uint32_t time_m);

    double error();
};

#endif