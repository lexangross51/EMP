#pragma once
#ifndef BASIS_FUNCTION_HPP
#define BASIS_FUNCTION_HPP

#include "area.hpp"
#include "head.hpp"

struct basis_function
{
    basis_function(omega& area)
    {
        _rk = area.start_point.r;
        _rk1 = area.end_point.r;
        _r_aver = (_rk + _rk1) / 2.0;

        _zk = area.start_point.z;
        _zk1 = area.end_point.z;
        _z_aver = (_zk + _zk1) / 2.0;

        _hr = _rk1 - _rk;
        _hz = _zk1 - _zk;
    }

    double psi(const uint8_t ifunc, const point2D& point)
    {
        double r = point.r;
        double z = point.z;

        switch (ifunc)
        {
        case 0:
            return 2.0 / (_hr * _hr) * (r - _r_aver) * (r - _rk1) * 2.0 / (_hz * _hz) * (z - _z_aver) * (z - _zk1);

        case 1:
            return (-4.0) / (_hr * _hr) * (r - _rk) * (r - _rk1) * 2.0 / (_hz * _hz) * (z - _z_aver) * (z - _zk1);

        case 2:
            return 2.0 / (_hr * _hr) * (r - _rk) * (r - _r_aver) * 2.0 / (_hz * _hz) * (z - _z_aver) * (z - _zk1);

        case 3:
            return 2.0 / (_hr * _hr) * (r - _r_aver) * (r - _rk1) * (-4.0) / (_hz * _hz) * (z - _zk) * (z - _zk1);

        case 4:
            return (-4.0) / (_hr * _hr) * (r - _rk) * (r - _rk1) * (-4.0) / (_hz * _hz) * (z - _zk) * (z - _zk1);

        case 5:
            return 2.0 / (_hr * _hr) * (r - _rk) * (r - _r_aver) * (-4.0) / (_hz * _hz) * (z - _zk) * (z - _zk1);

        case 6:
            return 2.0 / (_hr * _hr) * (r - _r_aver) * (r - _rk1) * 2.0 / (_hz * _hz) * (z - _zk) * (z - _z_aver);

        case 7:
            return (-4.0) / (_hr * _hr) * (r - _rk) * (r - _rk1) * 2.0 / (_hz * _hz) * (z - _zk) * (z - _z_aver);

        case 8:
            return 2.0 / (_hr * _hr) * (r - _rk) * (r - _r_aver) * 2.0 / (_hz * _hz) * (z - _zk) * (z - _z_aver);
        }
    }

    double d_psi(const uint8_t ifunc, const uint8_t ivar, const point2D& point)
    {
        double r = point.r;
        double z = point.z;

        switch (ivar)
        {
        case 1: {
            point2D next = point, prev = point;
            next.r += _hr;
            prev.r -= _hr;

            return (psi(ifunc, next) - psi(ifunc, prev)) / (2.0 * _hr);
        }

        case 2: {
            point2D next = point, prev = point;
            next.z += _hz;
            prev.z -= _hz;

            return (psi(ifunc, next) - psi(ifunc, prev)) / (2.0 * _hz);
        }
        }
    }

private:
    double _rk, _r_aver, _rk1;
    double _zk, _z_aver, _zk1;
    double _hr, _hz;
};

#endif