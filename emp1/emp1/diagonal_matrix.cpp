#include "diagonal_matrix.h"

diagonal_matrix::diagonal_matrix(const uint32_t _dim, const uint32_t _diags_count)
{
	dim = _dim;

	diags.resize(_diags_count);

	offset.resize(_diags_count);

	// ����� ���������� ��� ������� (������� �� � ��� ���)
	uint32_t under_main = (_diags_count - 1) / 2;

	offset[0] = 0;
	for (uint32_t i = 1; i < under_main + 1; i++)
	{
		offset[i] = -1 * i;
		offset[_diags_count - under_main + i - 1] = i;
	}

	for (uint32_t i = 0; i < _diags_count; i++)
	{
		if (i == 0)
			diags[0].resize(dim);
		else
			diags[i].resize(dim - abs(offset[i]));
	}
}

// ��������� ������� � ������� ������
void diagonal_matrix::to_dense(const std::string dir)
{
	// ���������� ���������� � ������(�������) ������������
	uint32_t under_main = (diags.size() - 1) / 2;

	matrix.resize(dim);
	for (int i = 0; i < dim; i++)
		matrix[i].resize(dim);

	// �������� �� ������� ������������, ���������� ������� ���������
	// �������� �������� ������������ ����� �������� ������������ 
	// ���������� ��������
	for (uint32_t i = 0; i < under_main + 1; i++)
	{
		// k - �������� �������� ��������� ������������ �������
		uint32_t k = abs(offset[i]);
		for (uint32_t j = 0; j < diags[i].size(); j++, k++)
		{
			matrix[k][j] = diags[i][j];

			if (i == 0)
				continue;

			matrix[j][k] = diags[under_main + i][j];
		}
	}

	//===============================================================
	std::ofstream dense(dir + "dense.txt");
	
	if (dense.is_open())
	{
		for (uint32_t i = 0; i < dim; i++)
		{
			for (uint32_t j = 0; j < dim; j++)
				dense << std::setw(5) << std::left << matrix[i][j];
			dense << std::endl << std::endl;
		}
		dense.close();
	}
	else
		throw "Can't open file";
}

// �������� ������� �� ������
std::unique_ptr<std::vector<double>> diagonal_matrix::dot(const std::vector<double>& vector)
{
	auto res = std::make_unique<std::vector<double>>(vector.size());

	// ����� ���������� � ������ (�������) ������������
	uint32_t under_main = (diags.size() - 1) / 2;

	for (uint32_t i = 0; i < under_main + 1; i++)
	{
		uint32_t k = abs(offset[i]);
		for (uint32_t j = 0; j < diags[i].size(); j++, k++)
		{
			(*res)[k] += diags[i][j] * vector[j];

			if (i == 0) continue;

			(*res)[j] += diags[under_main + i][j] * vector[k];
		}
	}

	return res;
}

// �������� ������ ������� �� ������
double diagonal_matrix::dot(const uint32_t row, const std::vector<double>& vector)
{
	double sum = 0.0;

	// ����� ���������� ���(���) �������
	uint32_t out_main = (diags.size() - 1) / 2;

	// ����� ���������� ��� �������, ������� ������ � row-�� ������
	uint32_t under_main = row > out_main ? out_main : row;

	// ����� ���������� ��� �������, ������� ������ � row-�� ������
	uint32_t above_main = dim - 1 - row > out_main ? out_main : dim - 1 - row;

	// �������� ������� ���������
	sum += diags[0][row] * vector[row];

	// �������� �������� ������ ����������
	for (uint32_t i = 1; i < under_main + 1; i++)
	{
		uint32_t k = abs(offset[i]);
		sum += diags[i][row - k] * vector[row - k];
	}

	// �������� �������� ������� ����������
	for (uint32_t i = 1; i < above_main + 1; i++)
	{
		uint32_t k = abs(offset[i]);
		sum += diags[out_main + k][row] * vector[row + i];
	}

	return sum;
}

// ������� ������������ ������������
void diagonal_matrix::make_dom()
{
	for (uint32_t row = 0; row < dim; row++)
	{
		diags[0][row] =  row == 0 ? 1.0 : 0.0;

		uint32_t out_main = (diags.size() - 1) / 2;
		uint32_t under_main = row > out_main ? out_main : row;
		uint32_t above_main = dim - 1 - row > out_main ? out_main : dim - 1 - row;

		for (uint32_t i = 1; i < under_main + 1; i++)
		{
			uint32_t k = abs(offset[i]);
			diags[0][row] += abs(diags[i][row - k]);
		}

		for (uint32_t i = 1; i < above_main + 1; i++)
		{
			uint32_t k = abs(offset[i]);
			diags[0][row] += abs(diags[out_main + k][row]);
		}
	}
}