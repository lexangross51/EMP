#include <iostream>
#include "mesh_generator.h"
#include "diagonal_matrix.h"
#include "solver.h"
#include "fdm.h"

int main()
{
	#pragma region Тесты с сеткой
	//mesh_generator mg;
	//std::vector<node> mesh;
	//mg.build_mesh(mesh, mesh_generator::grid_type::UNIFORM);
	#pragma endregion

	#pragma region Тесты с матрицей
	//diagonal_matrix dm(21, 7, 1);

	//for (uint32_t i = 0; i < dm.diags.size(); i++)
	//	for (uint32_t j = 0; j < dm.diags[i].size(); j++)
	//		dm.diags[i][j] = rand() % 10;

	//dm.make_dom();
	//dm.to_dense(directory);

	//std::vector<double> v(dm.dim);
	//std::vector<double> vec;
	//for (uint32_t i = 0; i < v.size(); i++)
	//	v[i] = i + 1;

	//auto res = dm.dot(v);
	//for (const auto& i : *res)
	//	std::cout << i << "\t";
	//std::cout << std::endl;

	//for (uint32_t i = 0; i < v.size(); i++)
	//	std::cout << dm.dot(i, v) << std::endl;

	#pragma endregion

	#pragma region Тесты с решателем
	//solver slv(10000, 1e-13);

	//for (double w = 0; w < 2; w += 0.1)
	//{
	//	auto result = slv.solve(w, dm, *res, vec);
	//	std::cout	<< std::left 
	//				<< std::setw(10) << w 
	//				<< std::setw(10) << result.first 
	//				<< std::setw(20) << result.second 
	//				<< std::endl;

	//	for (const auto& i : vec)
	//		std::cout << i << "\t";
	//	std::cout << std::endl;
	//}

	#pragma endregion

	#pragma region Тесты с классом fdm
	//fdm fdm;
	//fdm.mesh_to_slae(mesh);

	#pragma endregion


	return 0;
}