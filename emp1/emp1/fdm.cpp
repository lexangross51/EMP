#include "fdm.h"

void fdm::mesh_to_slae(mesh& mesh, function2D u, function2D f)
{
	A = new diagonal_matrix(mesh.dim(), 5, mesh.get_height());
	b.resize(mesh.dim());
	q.resize(mesh.dim());

	// Количество элементов по оси X
	uint32_t kx = mesh.get_width() + 1;

	for (uint32_t i = 0; i < mesh.dim(); i++)
	{
		switch (mesh[i].type)
		{
		// Внутренний узел
		case node::node_type::INTERNAL:
		{
			auto node = mesh[i];

			double hx = mesh[i + 1].p.x - node.p.x;
			double hx_prev = node.p.x - mesh[i - 1].p.x;

			double hy = mesh[i + kx].p.y - node.p.y;
			double hy_prev = node.p.y - mesh[i - kx].p.y;

			A->diags[0][i]		= +2.0 / (hx_prev * hx) + 2.0 / (hy_prev * hy) + node.gamma;
			A->diags[1][i - 1]	= -2.0 * node.lambda / (hx_prev * (hx + hx_prev));
			A->diags[3][i]		= -2.0 * node.lambda / (hx * (hx + hx_prev));
			A->diags[2][i - kx]	= -2.0 * node.lambda / (hy_prev * (hy + hy_prev));
			A->diags[4][i]		= -2.0 * node.lambda / (hy * (hy + hy_prev));

			b[i] = node.gamma * f(node.p.x, node.p.y);

			break;
		}

		// Фиктивный узел
		case node::node_type::FICTITIOUS:
		{
			A->diags[0][i] = 1.0;
			b[i] = 0.0;

			break;
		}

		// Граничный узел
		case node::node_type::BORDER:
		{
			switch (mesh[i].bc)
			{
			// Первые краевые
			case border::bound_cond::DIRICHLET:
			{
				A->diags[0][i] = 1.0;
				b[i] = u(mesh[i].p.x, mesh[i].p.y);
				break;
			}

			// Вторые краевые(+)
			case border::bound_cond::P_NEUMANN:
			{

				break;
			}
			// Вторые краевые(-)
			case border::bound_cond::M_NEUMANN:
			{

				break;
			}

			// Третьи краевые(+)
			case border::bound_cond::P_NEWTON:
			{

				break;
			}
			// Третьи краевые(-)
			case border::bound_cond::M_NEWTON:
			{

				break;
			}

			// Краевые не заданы
			case border::bound_cond::NONE:
			{
				break;
			}
			}
		}
		}
	}

	A->to_dense(directory);
}

// Решить системы и вывести результат
void fdm::calculate()
{
	solver slv(10000, 1e-13);

	for (double w = 0; w < 2; w += 0.1)
	{
		auto result = slv.solve(w, *A, b, q);

		write(directory, w, result);
	}
}

// Вывести результат в файл
void fdm::write(const std::string dir, const double w, const std::pair<uint32_t, double> result)
{
	std::ofstream res(dir + "result.txt", std::ios::app);
	
	if (res.is_open())
	{
		res << std::left
			<< std::setw(12) << "w: " << w << std::endl
			<< std::setw(12) << "iterations: " << result.first << std::endl
			<< std::setw(12) << "residual: " << result.second << std::endl;

		for (const auto& it : q)
			res << std::setw(5) << it;

		res << std::endl << std::endl;
	}
	else
		std::cerr << "Can't open file" << std::endl;
}