#include "fdm.h"

function2D u = [](double x, double y)
{
	return x + y;
};

void fdm::mesh_to_slae(std::vector<node>& mesh)
{
	A = new diagonal_matrix(mesh.size(), 5, 1);
	b.resize(mesh.size());
	q.resize(mesh.size());

	for (uint32_t i = 0; i < mesh.size(); i++)
	{
		switch (mesh[i].type)
		{
			// Внутренний узел
		case node::node_type::INTERNAL:
		{


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
			case border::bound_cond::DIRICHLET:
			{
				A->diags[0][i] = 1.0;
				b[i] = u(mesh[i].p.x, mesh[i].p.y);
				break;
			}

			case border::bound_cond::NEUMANN:
			{

				break;
			}

			case border::bound_cond::NEWTON:
			{

				break;
			}

			case border::bound_cond::NONE:
			{
				break;
			}

			break;
			}


		default:
			break;
		}
		}


	}
}