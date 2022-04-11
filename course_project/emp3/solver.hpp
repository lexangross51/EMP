#pragma once
#ifndef SOLVER_HPP
#define SOLVER_HPP

#include "matrix.hpp"
#include "utilities.hpp"
#include "decomposer.hpp"

struct result
{
    double residual;
    uint32_t iters;
};

enum class iterative_method
{
    LOS,
    LOS_LU,
    BCGSTAB_LU
};

class solver
{
    friend class mfe;

public:
    void init(
        const uint32_t max_iter, 
        const double eps, 
        const dvector& init_approx,
        iterative_method method
    );

    result solve(
        matrix& A, 
        dvector& b, 
        dvector& x
    );

private:
    iterative_method _method;

    // Параметры для итерационного решателя
    uint32_t _max_iter;
    double _eps;

    dvector _initial_approx;

    // Итерационные методы
    result LOS(matrix& A, dvector& b, dvector& x);

    void LU_direct(matrix& A, const dvector& b, dvector& x);
    void LU_reverse(matrix& A, const dvector& b, dvector& x);

    result LOS_LU(matrix& A, dvector& b, dvector& x);

    result BCGSTAB_LU(matrix& A, dvector& b, dvector& x);
};

#endif