#pragma once
#ifndef MESH_H
#define MESH_H

#include "head.h"

// �������� �������
struct finite_elem
{
    uint32_t number;

    double x;
    double x_next;

    double sigma;

    finite_elem(uint32_t _number, double _x, double _x_next, double _sigma) :
        number(_number), x(_x), x_next(_x_next), sigma(_sigma) {};
};

// �����
class mesh
{
public:
    enum class mesh_type
    {
        UNIFORM,
        NONUNIFORM
    };

public:
    inline void add_elem(const finite_elem& elem) { elems.push_back(elem); }
    inline finite_elem& operator[] (const uint32_t ielem) { return elems[ielem]; }

    inline void set_time(const std::vector<double>& _time) { T = _time; };
    inline double get_time(const uint32_t moment) { return T[moment]; }

    inline uint32_t get_fe_count(void) { return elems.size(); }
    inline uint32_t get_time_size(void) { return T.size(); }

    inline mesh_type get_area_type(void) { return area_type; }
    inline mesh_type get_time_type(void) { return time_type; }

    void save(const std::string dir = directory);

private:
    std::vector<finite_elem> elems;				// �������� ��������
    std::vector<double> T;						// ������� �������

    mesh_type area_type = mesh_type::UNIFORM;	// ��� ����� �� ������������
    mesh_type time_type = mesh_type::UNIFORM;	// ��� ����� �� �������
};

#endif