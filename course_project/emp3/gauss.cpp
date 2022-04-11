#include "gauss.hpp"

gauss::gauss()
{
    _gauss_points = { 0.0, sqrt(3.0 / 5.0), -sqrt(3.0 / 5.0) };
    _gauss_weights = { 8.0 / 9.0, 5.0 / 9.0, 5.0 / 9.0 };
}

gauss::~gauss()
{
    _gauss_points.clear();
    _gauss_weights.clear();
}

double gauss::integrate1D(const function2D& f, const omega& area)
{
    double r = area.start_point.r;
    double zk = area.start_point.z;
    double zk1 = area.end_point.z;

    double hz = abs(zk1 - zk);

    double pi, qi;

    double sum = 0.0;

    for (uint8_t i = 0; i < 3; i++)
    {
        qi = _gauss_weights[i];
        pi = (zk + zk1 + _gauss_points[i] * hz) / 2.0;

        sum += qi * f(r, pi);
    }

    return sum * hz / 2.0;
}

double gauss::integrate2D(const function2D& f, const omega& area)
{
    double rk = area.start_point.r;
    double rk1 = area.end_point.r;
    double zk = area.start_point.z;
    double zk1 = area.end_point.z;

    double hr = abs(rk1 - rk);
    double hz = abs(zk1 - zk);

    double pi, pj, qi, qj;

    double sum = 0.0;

    for (uint8_t i = 0; i < 3; i++)
    {
        qi = _gauss_weights[i];
        pi = (rk + rk1 + _gauss_points[i] * hr) / 2.0;

        for (uint8_t j = 0; j < 3; j++)
        {
            qj = _gauss_weights[j];
            pj = (zk + zk1 + _gauss_points[j] * hz) / 2.0;

            sum += qi * qj * f(pi, pj);
        }
    }

    return sum * hr * hz / 4.0;
}