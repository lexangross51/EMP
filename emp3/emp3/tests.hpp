#pragma once
#ifndef TESTS
#define TESTS

#include "mesh.hpp"
#include "mesh_generator.hpp"
#include "mfe.hpp"
#include "solver.hpp"
#include "utilities.hpp"

class tests
{
public:
	void run();

private:
	void omega_tests();

	void lambda_tests();

	void sigma_tests();

	void hi_tests();

	void calc_exact(space_grid& grid);

	void print_diff(dvector& solution, std::ofstream& file, std::string method);

	dvector exact;
};

#endif