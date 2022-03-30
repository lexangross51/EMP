#include "mfe.hpp"

mfe::mfe(space_grid& mesh)
{
	grid = &mesh;

	w = 0.0;

	A = new matrix(2 * grid->get_nodes_count());

	local_p.resize(8);
	local_c.resize(8);

	local_vec_fs.resize(8);
	local_vec_fc.resize(8);

	for (uint8_t i = 0; i < 8; i++)
	{
		local_p[i].resize(8);
		local_c[i].resize(8);
	}
}

void mfe::build_local_matrices(finite_elem& elem)
{
	omega omega = { point3D(grid->get_point(elem[0])), point3D(grid->get_point(elem[7])) };

	basis_function bf(omega);

	function3D f;

	for (uint8_t i = 0; i < 8; i++)
	{
		for (uint8_t j = 0; j <= i; j++)
		{
			f = [&](double x, double y, double z)
			{
				point3D point(x, y, z);

				double psi_i = bf.psi(i, point);
				double psi_j = bf.psi(j, point);

				double d_psi_i_x = bf.d_psi(i, 1, point);
				double d_psi_j_x = bf.d_psi(j, 1, point);

				double d_psi_i_y = bf.d_psi(i, 2, point);
				double d_psi_j_y = bf.d_psi(j, 2, point);

				double d_psi_i_z = bf.d_psi(i, 3, point);
				double d_psi_j_z = bf.d_psi(j, 3, point);

				return (elem.lambda * (
					d_psi_i_x * d_psi_j_x + 
					d_psi_i_y * d_psi_j_y + 
					d_psi_i_z * d_psi_j_z
				) - w * w * elem.hi * psi_i * psi_j);
			};

			local_p[i][j] = local_p[j][i] = gauss.integrate3D(f, omega);


			// -----------------------------------------------------------
			f = [&](double x, double y, double z)
			{
				point3D point(x, y, z);

				double psi_i = bf.psi(i, point);
				double psi_j = bf.psi(j, point);

				return w * elem.sigma * psi_i * psi_j;
			};

			local_c[i][j] = local_c[j][i] = gauss.integrate3D(f, omega);
		}
	}
}

void mfe::build_local_vector(finite_elem& elem, r_function3D& fs, r_function3D& fc)
{
	omega omega = { point3D(grid->get_point(elem[0])), point3D(grid->get_point(elem[7])) };

	basis_function bf(omega);

	function3D f;

	for (uint8_t i = 0; i < 8; i++)
	{
		f = [&](double x, double y, double z)
		{
			point3D point(x, y, z);

			double psi_i = bf.psi(i, point);

			return fs(x, y, z, elem.lambda, elem.sigma, elem.hi, w) * psi_i;
		};

		local_vec_fs[i] = gauss.integrate3D(f, omega);


		// -----------------------------------------------------------
		f = [&](double x, double y, double z)
		{
			point3D point(x, y, z);

			double psi_i = bf.psi(i, point);

			return fc(x, y, z, elem.lambda, elem.sigma, elem.hi, w) * psi_i;
		};

		local_vec_fc[i] = gauss.integrate3D(f, omega);
	}
}

void mfe::add_to_global_matrix(const uint32_t i, const uint32_t j, const double val)
{
	if (i == j)
	{
		A->di[i] += val;
		return;
	}

	if (i < j)
	{
		for (uint32_t ind = A->ig[j]; ind < A->ig[j + 1]; ind++)
		{
			if (A->jg[ind] == i)
			{
				A->ggu[ind] += val;
				return;
			}
		}
	}
	else
	{
		for (uint32_t ind = A->ig[i]; ind < A->ig[i + 1]; ind++)
		{
			if (A->jg[ind] == j)
			{
				A->ggl[ind] += val;

			}
		}
	}
}

void mfe::assembly_global_matrix_and_vector(r_function3D& fs, r_function3D& fc)
{
	portrait_generator::portrait(*grid, A->ig, A->jg);

	A->di.resize(A->dim);
	A->ggl.resize(A->ig.back());
	A->ggu.resize(A->ig.back());
	b.resize(A->dim);

	for (uint32_t ielem = 0; ielem < grid->get_elems_count(); ielem++)
	{
		finite_elem elem = grid->get_elem(ielem);

		build_local_matrices(elem);

		#pragma region Проверка
		//std::cout.precision(3);

		//std::cout << "LOCAL P\n";
		//for (uint32_t j = 0; j < 8; j++)
		//{ 
		//	for (uint32_t k = 0; k < 8; k++)
		//		std::cout << local_p[j][k] << "    ";
		//	std::cout << std::endl;
		//}

		//std::cout << "LOCAL C\n";
		//for (uint32_t j = 0; j < 8; j++)
		//{
		//	for (uint32_t k = 0; k < 8; k++)
		//		std::cout << local_c[j][k] << "    ";
		//	std::cout << std::endl;
		//}
		//std::cout << "---------------------------------------------------------------------\n";
		//std::cout << std::endl;
		#pragma endregion

		// ------------------------------------------------------------
		build_local_vector(elem, fs, fc);

		#pragma region Проверка
		//std::cout.precision(3);

		//std::cout << "LOCAL FS\n";
		//for (uint32_t j = 0; j < 8; j++)
		//{ 
		//	std::cout << local_vec_fs[j] << "    ";
		//}

		//std::cout << std::endl;
		//std::cout << "LOCAL FC\n";
		//for (uint32_t j = 0; j < 8; j++)
		//{
		//	std::cout << local_vec_fc[j] << "    ";
		//}
		//std::cout << std::endl;
		//std::cout << "---------------------------------------------------------------------\n";
		//std::cout << std::endl;
		#pragma endregion

		for (uint8_t i = 0; i < 8; i++)
		{
			b[2 * elem[i]] += local_vec_fs[i];
			b[2 * elem[i] + 1] += local_vec_fc[i];

			for (uint8_t j = 0; j < 8; j++)
			{
				add_to_global_matrix(2 * elem[i], 2 * elem[j], local_p[i][j]);
				add_to_global_matrix(2 * elem[i] + 1, 2 * elem[j] + 1, local_p[i][j]);
				add_to_global_matrix(2 * elem[i] + 1, 2 * elem[j], local_c[i][j]);
				add_to_global_matrix(2 * elem[i], 2 * elem[j] + 1, -local_c[i][j]);
			}
		}
	}
}