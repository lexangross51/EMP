#include <iostream>
#include "mesh_generator.h"
#include "mesh.h"
#include "diagonal_matrix.h"
#include "solver.h"
#include "fdm.h"

func2D u = [](double x, double y)
{
	return x + y;
};

function2D f = [](double x, double y, double coef)
{
	return coef * (x + y);
};

int main()
{
	mesh_generator mg;
	mesh mesh;

	mg.build_mesh(mesh, mesh::mesh_type::NONUNIFORM);

	fdm fdm;
	fdm.mesh_to_slae(mesh, u, f);
	fdm.calculate();

	return 0;
}