#include <iostream>
#include "mesh_generator.h"

int main()
{
	mesh_generator mg;
	std::vector<node> mesh;

	mg.build_mesh(mesh);

	return 0;
}