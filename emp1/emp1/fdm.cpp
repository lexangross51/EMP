#include "fdm.h"

const double RELAXATION_PARAMETER = 1.26;

// Преобразовываем сетку в СЛАУ
void fdm::mesh_to_slae(mesh& mesh, func2D_u& u, func2D_f& f)
{
	A = new diagonal_matrix(mesh.dim(), 5, mesh.get_width() - 1);
	b.resize(mesh.dim());

	// Количество элементов по оси X
	uint32_t kx = mesh.get_width() + 1;

	for (uint32_t i = 0; i < mesh.dim(); i++)
	{
		auto node = mesh[i];

		// Внутренний узел
		if (node.type == 0)
		{
			double hx = mesh[i + 1].p.x - node.p.x;
			double hx_prev = node.p.x - mesh[i - 1].p.x;

			double hy = mesh[i + kx].p.y - node.p.y;
			double hy_prev = node.p.y - mesh[i - kx].p.y;

			A->diags[0][i]		= (+ 2.0 / (hx_prev * hx) + 2.0 / (hy_prev * hy)) * node.lambda + node.gamma;
			A->diags[1][i - 1]	= -2.0 * node.lambda / (hx_prev * (hx + hx_prev));
			A->diags[3][i]		= -2.0 * node.lambda / (hx * (hx + hx_prev));
			A->diags[2][i - kx]	= -2.0 * node.lambda / (hy_prev * (hy + hy_prev));
			A->diags[4][i]		= -2.0 * node.lambda / (hy * (hy + hy_prev));

			b[i] = f(node.p.x, node.p.y, node.lambda, node.gamma);
		}

		// Фиктивный узел
		else if (node.type == -1)
		{
			A->diags[0][i] = 1.0;
			b[i] = 0.0;
		}

		// Граничный узел
		else if (node.type == 1 || node.type == 2 || node.type == 3 || node.type == 4)
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

			// Вторые краевые
			case border::bound_cond::NEUMANN:
			{
				double hx, hy;
				if (node.type == 2)
				{
					hx = mesh[i + 1].p.x - node.p.x;
					A->diags[0][i] = -1.0 / hx;
					A->diags[3][i] = 1.0 / hx;
					b[i] = node.lambda * du_dx(u, node.p);
				}
				else if (node.type == 4)
				{
					hx = node.p.x - mesh[i - 1].p.x;
					A->diags[0][i] = 1.0 / hx;
					A->diags[1][i - 1] = -1.0 / hx;
					b[i] = node.lambda * du_dx(u, node.p);
				}
				else if (node.type == 1)
				{
					hy = mesh[i + kx].p.y - node.p.y;
					A->diags[0][i] = -1.0 / hy;
					A->diags[4][i] = 1.0 / hy;
					b[i] = node.lambda * du_dy(u, node.p);
				}
				else if (node.type == 3)
				{
					hy = node.p.y - mesh[i - kx].p.y;
					A->diags[0][i] = 1.0 / hy;
					A->diags[2][i - kx] = -1.0 / hy;
					b[i] = node.lambda * du_dy(u, node.p);
				}
				break;
			}

			// Третьи краевые
			case border::bound_cond::NEWTON:
			{
				double hx, hy;
				if (node.type == 2)
				{
					hx = mesh[i + 1].p.x - node.p.x;
					A->diags[0][i] = -1.0 / hx + beta;
					A->diags[3][i] = 1.0 / hx;
					b[i] = node.lambda * du_dx(u, node.p) + beta * (u(node.p.x, node.p.y) - u_beta);
				}
				else if (node.type == 4)
				{
					hx = node.p.x - mesh[i - 1].p.x;
					A->diags[0][i] = 1.0 / hx + beta;
					A->diags[1][i - 1] = -1.0 / hx;
					b[i] = node.lambda * du_dx(u, node.p) + beta * (u(node.p.x, node.p.y) - u_beta);
				}
				else if (node.type == 1)
				{
					hy = mesh[i + kx].p.y - node.p.y;
					A->diags[0][i] = -1.0 / hy + beta;
					A->diags[4][i] = 1.0 / hy;
					b[i] = node.lambda * du_dy(u, node.p) + beta * (u(node.p.x, node.p.y) - u_beta);
				}
				else if (node.type == 3)
				{
					hy = node.p.y - mesh[i - kx].p.y;
					A->diags[0][i] = 1.0 / hy + beta;
					A->diags[2][i - kx] = -1.0 / hy;
					b[i] = node.lambda * du_dy(u, node.p) + beta * (u(node.p.x, node.p.y) - u_beta);
				}
				break;
			}

			// По-умолчанию задаем первые краевые
			default:
			{
				A->diags[0][i] = 1.0;
				b[i] = u(mesh[i].p.x, mesh[i].p.y);
				break;
			}
			}
		}
	}
}

// Решить систему и вернуть результат
std::pair<uint32_t, double> fdm::calculate(std::vector<double>& q)
{
	q.resize(b.size());

	solver slv(100000, 1e-14);

	return slv.solve(RELAXATION_PARAMETER, *A, b, q);
}

// 1-я производная по x
double fdm::du_dx(func2D_u& f, point& p)
{
	double h = 1e-5;
	return (f(p.x + h, p.y) - f(p.x, p.y)) / h;
}

// 1-я производная по y
double fdm::du_dy(func2D_u& f, point& p)
{
	double h = 1e-5;
	return (f(p.x, p.y + h) - f(p.x, p.y)) / h;
}