#include "mesh.hpp"
#include "mesh_generator.hpp"
#include "mfe.hpp"

function3D us = [](double x, double y, double z) { return x*x*x + y*y*y + z*z*z; };
function3D uc = [](double x, double y, double z) { return 2*x*x*x - y*y*y + 3*z*z*z; };

r_function3D fs = [](double x, double y, double z, double lambda, double sigma, double hi, double w)
{
	return -lambda * 6 * (x + y + z) - w * sigma * uc(x, y, z) - w * w * hi * us(x, y, z);
};

r_function3D fc = [](double x, double y, double z, double lambda, double sigma, double hi, double w)
{
	return -lambda * 6 * (2 * x - y + 3 * z) + w * sigma * us(x, y, z) - w * w * hi * uc(x, y, z);
};

int main()
{
	space_grid_generator sgg;
	space_grid* sg = nullptr;

	sgg.build_mesh(sg);

	mfe mfe(*sg);


	mfe.set_w(2);
	mfe.assembly_global_matrix_and_vector(fs, fc);


	return 0;
}