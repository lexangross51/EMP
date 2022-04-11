#pragma once
#ifndef MESH_GENERATOR_HPP
#define MESH_GENERATOR_HPP

#include "mesh.hpp"
#include "area.hpp"

// ---------------------------------------------------------------------------
// Абстрактный класс генератора сетки
class grid_generator
{
protected:
    virtual void read_data(std::string path = directory + "input_data\\") = 0;

    virtual void generate_nodes() = 0;
};


// =================== ГЕНЕРАЦИЯ СЕТКИ ПО ПРОСТРАНСТВУ ====================
class space_grid_generator : public grid_generator
{
public:
    space_grid_generator();

    void build_mesh(space_grid*& grid);

public:
    // Входные данные для построения сетки
    double _rw, _rb;
    double _depth;
    dvector _heights;
    dvector _lambda_s, _sigma_s, _hi_s;

    uint32_t _nr;
    std::vector<uint32_t> _nz;

    double _kr;
    std::vector<double> _kz;

    mesh_type _type;
    uint8_t _nested;

    std::array<boundary_cond, 4> _bconds;

    // -----------------------------------
    std::set<double> _r, _z;

    void read_data(std::string path = directory + "input_data\\") override;

    void generate_nodes() override;

    void make_bc_nodes(space_grid& grid);
};


// ====================== ГЕНЕРАЦИЯ СЕТКИ ПО ВРЕМЕНИ =======================
class time_grid_generator : public grid_generator
{
public:
    time_grid_generator();

    void build_mesh(time_grid*& grid);

private:
    // -------------------
    dvector _time_layers;
    
    double _time_begin, _time_end;
    double _nt;
    double _kt;
    mesh_type _type;
    uint8_t _nested;

    // -------------------
    void read_data(std::string path = directory + "input_data\\") override;

    void generate_nodes() override;
};


#endif