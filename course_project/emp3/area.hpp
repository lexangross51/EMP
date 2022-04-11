#pragma once 
#ifndef AREA_HPP
#define AREA_HPP

struct point2D
{
    double r, z;

    point2D(double _r = 0.0, double _z = 0.0) :
        r(_r), z(_z) {};
};

// Представляет собой границы конечного элемента
struct omega
{
    point2D start_point;
    point2D end_point;
};

#endif