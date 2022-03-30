#pragma once
#ifndef BASIS_FUNCTION_HPP
#define BASIS_FUNCTION_HPP

#include "finite_element.hpp"
#include "area.hpp"

struct basis_function
{
	// Границы области, в которой определены базисные функции
	double xk, xk1;
	double yk, yk1;
	double zk, zk1;
	double hx, hy, hz;

	basis_function(omega& _area) {
		xk = _area.start_point.x;
		yk = _area.start_point.y;
		zk = _area.start_point.z;
		xk1 = _area.end_point.x;
		yk1 = _area.end_point.y;
		zk1 = _area.end_point.z;

		hx = xk1 - xk;
		hy = yk1 - yk;
		hz = zk1 - zk;
	}

	double psi(uint32_t ifunc, point3D& point)
	{
		double x = point.x;
		double y = point.y;
		double z = point.z;

		switch (ifunc)
		{
		case 0:
			return (xk1 - x) * (yk1 - y) * (zk1 - z) / (hx * hy * hz);

		case 1:
			return (x - xk) * (yk1 - y) * (zk1 - z) / (hx * hy * hz);

		case 2:
			return (xk1 - x) * (y - yk) * (zk1 - z) / (hx * hy * hz);

		case 3:
			return (x - xk) * (y - yk) * (zk1 - z) / (hx * hy * hz);

		case 4:
			return (xk1 - x) * (yk1 - y) * (z - zk) / (hx * hy * hz);

		case 5:
			return (x - xk) * (yk1 - y) * (z - zk) / (hx * hy * hz);

		case 6:
			return (xk1 - x) * (y - yk) * (z - zk) / (hx * hy * hz);

		case 7:
			return (x - xk) * (y - yk) * (z - zk) / (hx * hy * hz);
		}
	}

	double d_psi(uint8_t ifunc, uint8_t ivar, point3D& point)
	{
		double h = 1e-7;

		switch (ivar)
		{
		case 1: {
			point3D next = point, prev = point;
			next.x += h;
			prev.x -= h;
			return (psi(ifunc, next) - psi(ifunc, prev)) / (2 * h);
		}
		case 2: {
			point3D next = point, prev = point;
			next.y += h;
			prev.y -= h;
			return (psi(ifunc, next) - psi(ifunc, prev)) / (2 * h);
		}
		case 3: {
			point3D next = point, prev = point;
			next.z += h;
			prev.z -= h;
			return (psi(ifunc, next) - psi(ifunc, prev)) / (2 * h);
		}
		}
	}
};

#endif