#include "matrix.hpp"

std::unique_ptr<dvector> matrix::dot(const dvector& vector)
{
    auto result = std::make_unique <std::vector <double>>(vector.size());

    for (uint32_t i = 0; i < result->size(); i++)
    {
        (*result)[i] = _di[i] * vector[i];
        for (uint32_t j = _ig[i]; j < _ig[i + 1]; j++)
        {
            (*result)[i] += _ggl[j] * vector[_jg[j]];
            (*result)[_jg[j]] += _ggu[j] * vector[i];
        }
    }

    return result;
}

void matrix::to_dense()
{
    _dense_matrix.resize(_dim);

    for (uint32_t i = 0; i < _dense_matrix.size(); i++)
        _dense_matrix[i].resize(_dim, 0);

    for (uint32_t i = 0; i < _dim; i++)
    {
        _dense_matrix[i][i] = _di[i];
        for (uint32_t j = _ig[i]; j < _ig[i + 1]; j++)
        {
            _dense_matrix[i][_jg[j]] = _ggl[j];
            _dense_matrix[_jg[j]][i] = _ggu[j];
        }
    }
}

void matrix::save(matrix_type type, std::string path)
{
    switch (type)
    {
    case matrix_type::DENSE: {
        std::ofstream dense_out(path + "dense.txt");

        dense_out.precision(3);

        for (uint32_t i = 0; i < _dim; i++)
        {
            for (uint32_t j = 0; j < _dim; j++)
                dense_out << std::setw(10) << _dense_matrix[i][j];
            dense_out << std::endl << std::endl;
        }

        dense_out.close();
    }
    break;

    case matrix_type::SPARSE: {
        std::ofstream sparse_out(path + "sparse.json");
        nlohmann::json sparse{};

        sparse["ig"] = _ig;
        sparse["jg"] = _jg;
        sparse["di"] = _di;
        sparse["ggl"] = _ggl;
        sparse["ggu"] = _ggu;

        sparse_out << sparse;
        sparse_out.close();
    }
    break;

    default:
        throw "What type is this?\n";
    }
}

void matrix::null()
{
    for (uint32_t i = 0; i < _di.size(); i++)
        _di[i] = 0.0;

    for (uint32_t i = 0; i < _ggl.size(); i++)
    {
        _ggl[i] = 0.0;
        _ggu[i] = 0.0;
    }
}