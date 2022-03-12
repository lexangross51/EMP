#pragma once
#ifndef MESH_GENERATOR_H
#define MESH_GENERATOR_H

#include "head.h"
#include "mesh.h"

class mesh_generator
{
public:
    mesh_generator()
    {
        area_start = area_end = 0; 
        nx = kx = 0;
        area_nested = 0;

        time_start = time_end = 0; 
        nt = kt = 0;
        time_nested = 0;

        sigma = 0.0;
    }

    void build_mesh(mesh& mesh, mesh::mesh_type area_type, mesh::mesh_type time_type);

private:
    double area_start;
    double area_end;

    double sigma;

    uint32_t nx;
    double kx;

    uint32_t area_nested;

    std::vector<double> X;

    //----------------------
    double time_start;
    double time_end;

    uint32_t nt;
    double kt;

    uint32_t time_nested;

    std::vector<double> T;

    void input(const std::string dir = directory);

    void generate_XT(mesh::mesh_type area_type, mesh::mesh_type time_type);
};

#endif