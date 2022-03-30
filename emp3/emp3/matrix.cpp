#include "matrix.hpp"

std::unique_ptr<std::vector<double>> matrix::dot(const std::vector<double>& vector)
{
	auto result = std::make_unique <std::vector <double>>(vector.size());

	for (uint32_t i = 0; i < result->size(); i++)
	{
		(*result)[i] = di[i] * vector[i];
		for (uint32_t j = ig[i]; j < ig[i + 1]; j++)
		{
			(*result)[i] += ggl[j] * vector[jg[j]];
			(*result)[jg[j]] += ggu[j] * vector[i];
		}
	}

	return result;
}

void matrix::to_dense()
{
	dense_matrix.resize(dim);

	for (uint32_t i = 0; i < dense_matrix.size(); i++)
		dense_matrix[i].resize(dim, 0);

	for (uint32_t i = 0; i < dim; i++)
	{
		dense_matrix[i][i] = di[i];
		for (uint32_t j = ig[i]; j < ig[i + 1]; j++)
		{
			dense_matrix[i][jg[j]] = ggl[j];
			dense_matrix[jg[j]][i] = ggu[j];
		}
	}
}

void matrix::to_sparse()
{
	di.resize(dense_matrix.size(), 0);
	ig.resize(dense_matrix.size() + 1, 0);

	for (uint32_t i = 0; i < di.size(); i++)
		di[i] = dense_matrix[i][i];

	ig[0] = ig[1] = 0;
	for (uint32_t i = 1; i < dense_matrix.size(); i++)
	{
		uint32_t count = 0;
		for (uint32_t j = 0; j < i; j++)
		{
			if (dense_matrix[i][j] != 0.0 && dense_matrix[j][i] != 0.0 ||
				dense_matrix[i][j] == 0.0 && dense_matrix[j][i] != 0.0 ||
				dense_matrix[i][j] != 0.0 && dense_matrix[j][i] == 0.0)
				count++;
		}
		ig[i + 1] = ig[i] + count;
	}

	jg.resize(ig.back(), 0);
	uint32_t k = 0;
	for (uint32_t i = 0; i < dense_matrix.size(); i++)
	{
		for (uint32_t j = 0; j < i; j++)
		{
			if (dense_matrix[i][j] != 0.0 && dense_matrix[j][i] != 0.0 ||
				dense_matrix[i][j] == 0.0 && dense_matrix[j][i] != 0.0 ||
				dense_matrix[i][j] != 0.0 && dense_matrix[j][i] == 0.0)
			{
				jg[k] = j;
				k++;
			}
		}
	}

	ggl.resize(ig.back(), 0);
	ggu.resize(ig.back(), 0);
	k = 0;
	for (uint32_t i = 0; i < dim; i++)
	{
		for (uint32_t j = 0; j < i; j++)
		{
			if (dense_matrix[i][j] != 0.0 && dense_matrix[j][i] != 0.0 ||
				dense_matrix[i][j] == 0.0 && dense_matrix[j][i] != 0.0 ||
				dense_matrix[i][j] != 0.0 && dense_matrix[j][i] == 0.0)
			{
				ggl[k] = dense_matrix[i][j];
				ggu[k] = dense_matrix[j][i];
				k++;
			}
		}
	}
}

void matrix::to_profile()
{
	uint32_t count = 0;
	bool first = true;
	uint32_t ggl_size = 0;
	uint32_t k = 0;

	di.resize(dim);
	ig.resize(dim + 1);

	jg[0] = jg[1] = 0;

	for (uint32_t i = 0; i < dim; i++)
	{
		di[i] = dense_matrix[i][i];

		for (uint32_t j = 0; j < i; j++)
		{
			if (dense_matrix[i][j] == 0.0 && first && dense_matrix[j][i] != 0.0 ||
				dense_matrix[i][j] != 0.0 && first && dense_matrix[j][i] == 0.0 ||
				dense_matrix[i][j] != 0.0 && first && dense_matrix[j][i] != 0.0)
			{
				count++;
				first = false;
			}
			else if (!first)
				count++;
		}
		jg[i + 1] = jg[i] + count;
		count = 0;
		first = true;
	}

	ggl.resize(ig.back(), 0);
	ggu.resize(ig.back(), 0);
	for (uint32_t i = 0; i < dim; i++)
	{
		for (uint32_t j = 0; j < i; j++)
		{
			if (dense_matrix[i][j] == 0.0 && first && dense_matrix[j][i] != 0.0 ||
				dense_matrix[i][j] != 0.0 && first && dense_matrix[j][i] == 0.0 ||
				dense_matrix[i][j] != 0.0 && first && dense_matrix[j][i] != 0.0)
			{
				ggl[k] = dense_matrix[i][j];
				ggu[k] = dense_matrix[j][i];
				k++;

				first = false;
			}
			else if (!first)
			{
				ggl[k] = dense_matrix[i][j];
				ggu[k] = dense_matrix[j][i];
				k++;
			}
		}
		first = true;
	}
}

void matrix::save(matrix_type type, std::string path)
{
	switch (type)
	{
	case matrix_type::DENSE: {
		std::ofstream dense_out(path + "dense.json");
		nlohmann::json dense{};

		dense["dense"] = dense_matrix;

		dense_out << dense;
		dense_out.close();
	}
	break;

	case matrix_type::PROFILE: {
		std::ofstream profile_out(path + "profile.json");
		nlohmann::json profile{};

		profile["ig"] = ig;
		profile["di"] = di;
		profile["ggl"] = ggl;
		profile["ggu"] = ggu;

		profile_out << profile;
		profile_out.close();
	}
	break;

	case matrix_type::SPARSE: {
		std::ofstream sparse_out(path + "sparse.json");
		nlohmann::json sparse{};

		sparse["ig"] = ig;
		sparse["jg"] = jg;
		sparse["di"] = di;
		sparse["ggl"] = ggl;
		sparse["ggu"] = ggu;

		sparse_out << sparse;
		sparse_out.close();
	}
	break;

	default:
		throw "What type is this?\n";
	}
}