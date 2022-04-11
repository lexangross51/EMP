#pragma once
#ifndef GAUSS_HPP
#define GAUSS_HPP

#include "head.hpp"
#include "area.hpp"

class gauss
{
public:
    gauss();
    ~gauss();

    double integrate1D(const function2D& f, const omega& area);

    double integrate2D(const function2D& f, const omega& area);

private:
    dvector _gauss_points;
    dvector _gauss_weights;
};

#endif