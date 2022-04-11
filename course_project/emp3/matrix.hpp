#pragma once
#ifndef MATRIX_HPP
#define MATRIX_HPP

#include "head.hpp"

enum class matrix_type
{
    DENSE,
    PROFILE,
    SPARSE
};

class matrix
{
    friend class solver;
    friend class decomposer;
    friend class mfe;

public:
    matrix(uint32_t _size) 
    { 
        _dim = _size; 
        _di.resize(_dim);
    }

    std::unique_ptr<dvector> dot(const dvector& vector);

    inline uint32_t size() const { return _dim; };

    inline uint32_t ig(const uint32_t index) const { return _ig[index]; }
    inline uint32_t jg(const uint32_t index) const { return _jg[index]; }
    inline double di(const uint32_t index) const { return _di[index]; }
    inline double ggl(const uint32_t index) const { return _ggl[index]; }
    inline double ggu(const uint32_t index) const { return _ggu[index]; }

    inline void di_s(const uint32_t index, const double value) { _di[index] = value; }
    inline void ggl_s(const uint32_t index, const double value) { _ggl[index] = value; }
    inline void ggu_s(const uint32_t index, const double value) { _ggu[index] = value; }

    inline void di_a(const uint32_t index, const double value) { _di[index] += value; }
    inline void ggl_a(const uint32_t index, const double value) { _ggl[index] += value; }
    inline void ggu_a(const uint32_t index, const double value) { _ggu[index] += value; }

    inline void set_portrait(const std::vector<uint32_t>& ig, const std::vector<uint32_t>& jg)
    {
        _ig = ig;
        _jg = jg;
        _ggl.resize(_ig.back());
        _ggu.resize(_ig.back());
    }

    void null();

    void to_dense();
    void save(matrix_type type, std::string path = directory + "matrix\\");

private:
    uint32_t _dim;

    std::vector<uint32_t> _ig, _jg;
    dvector _di, _ggl, _ggu;

    std::vector<dvector> _dense_matrix;
};

#endif