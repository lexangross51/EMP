#include "mesh.h"

void mesh::save(std::string path)
{
	std::ofstream grid_out("mesh.json");
	nlohmann::json grid{};

	for (uint32_t i = 0; i < dimension; i++)
	{

	}


}