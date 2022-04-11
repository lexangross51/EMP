#pragma once
#ifndef MESH_HPP
#define MESH_HPP

#include "finite_element.hpp"
#include "boundary.hpp"
#include "area.hpp"

enum class mesh_type
{
    UNIFORM,
    NONUNIFORM
};


// ==================== ¿¡—“–¿ “Õ€…  À¿—— —≈“ » =====================
class grid
{
public:
    grid() : _type(mesh_type::UNIFORM) {};

    virtual void set_type(const mesh_type type) { _type = type; }

    virtual inline uint32_t get_nodes_count(void) const = 0;
    virtual inline mesh_type get_type(void) const { return _type; }

    virtual void save(std::string path = directory + "mesh\\") = 0;

protected:
    mesh_type _type;
};


// ===================== —≈“ ¿ œŒ œ–Œ—“–¿Õ—“¬” ======================
class space_grid : public grid
{
    friend class space_grid_generator;

public:
    space_grid(const uint32_t elems_cnt = 0, const uint32_t points_cnt = 0) {
        _nr = _nz = 0;
        _elems.resize(elems_cnt);
        _points.resize(points_cnt);
    };

    ~space_grid() {
        _elems.clear();
        _points.clear();
        _dirichlet.clear();
        _neumann.clear();
    }

    inline void set_width(const uint32_t nr) { _nr = nr; }
    inline void set_height(const uint32_t nz) { _nz = nz; }

    inline void add_point(const point2D& point, const uint32_t pos) { _points[pos] = point; }
    inline void add_elem(const finite_element& element, const uint32_t pos) { _elems[pos] = element; }

    inline void set_dirichlet_conds(std::vector<dirichlet_cond>& nodes) { _dirichlet = nodes; }
    inline void set_neumann_conds(std::vector<neumann_cond>& nodes) { _neumann = nodes; }

    inline uint32_t get_width(void) const { return _nr; }
    inline uint32_t get_height(void) const { return _nz; }

    inline uint32_t get_elems_count(void) const { return _elems.size(); }
    inline uint32_t get_nodes_count(void) const override { return _points.size(); }

    inline finite_element& get_elem(const uint32_t index) { return _elems[index]; }
    inline point2D& get_point(const uint32_t index) { return _points[index]; }

    inline uint32_t get_first_bound_count(void) const { return _dirichlet.size(); }
    inline uint32_t get_second_bound_count(void) const { return _neumann.size(); }

    inline dirichlet_cond& get_dirichlet_cond(const uint32_t num) { return _dirichlet[num]; }
    inline neumann_cond& get_neumann_cond(const uint32_t num) { return _neumann[num]; }

    void save(std::string path = directory + "mesh\\") override;

private:
    uint32_t _nr;
    uint32_t _nz;

    std::vector<finite_element> _elems;
    std::vector<point2D> _points;

    std::vector<dirichlet_cond> _dirichlet;
    std::vector<neumann_cond> _neumann;
};


// ======================== —≈“ ¿ œŒ ¬–≈Ã≈Õ» =========================
class time_grid : public grid
{
    friend class time_grid_generator;

public:
    time_grid(const uint32_t time_layers = 0) {
        _layers.resize(time_layers);
    }

    ~time_grid() {
        _layers.clear();
    }

    inline double get_time(const uint32_t layer_pos) { return _layers[layer_pos]; }
    inline uint32_t get_nodes_count(void) const override { return _layers.size(); }

    void save(std::string path = directory + "mesh\\") override;

private:
    dvector _layers;
};

#endif