#include "solver.hpp"

void solver::init(const uint32_t max_iter, const double eps, const dvector& init_approx, iterative_method method)
{
    _max_iter = max_iter;
    _eps = eps;
    _initial_approx = init_approx;

    for (uint32_t i = 0; i < _initial_approx.size(); i++)
        _initial_approx[i] = 1;

    _method = method;
}

result solver::solve(matrix& A, dvector& b, dvector& x)
{
    if (_method == iterative_method::LOS)
        return LOS(A, b, x);

    if (_method == iterative_method::LOS_LU)
        return LOS_LU(A, b, x);

    if (_method == iterative_method::BCGSTAB_LU)
        return BCGSTAB_LU(A, b, x);
}

void solver::LU_direct(matrix& A, const dvector& b, dvector& x)
{
    x = b;

    for (uint32_t i = 0; i < x.size(); i++)
    {
        double sum = 0.0;

        for (uint32_t j = A._ig[i]; j < A._ig[i + 1]; j++)
            sum += A._ggl[j] * x[A._jg[j]];

        x[i] -= sum;
        x[i] /= A._di[i];
    }
}

void solver::LU_reverse(matrix& A, const dvector& b, dvector& x)
{
    x = b;

    for (__int64 i = A.size() - 1; i >= 0; i--)
    {
        for (uint32_t j = A._ig[i]; j < A._ig[i + 1]; j++)
            x[A._jg[j]] -= A._ggu[j] * x[i];
    }
}

result solver::LOS_LU(matrix& _A, dvector& b, dvector& x)
{
    matrix A = _A;

    result res;

    uint32_t dim = _A.size();

    dvector r(dim), z(dim), p(dim);
    dvector LAU(dim), U(dim);

    double alpha, beta;

    x = _initial_approx;

    decomposer::LU(A);

    LU_direct(A, b - *_A.dot(x), r);
    LU_reverse(A, r, z);
    LU_direct(A, *_A.dot(z), p);


    uint32_t k = 0;
    double rr = (r * r);

    for (; k < _max_iter && rr >= _eps; k++)
    {
        double pp = p * p;
        alpha = (p * r) / pp;

        x = x + alpha * z;

        rr = (r * r);

        r = r - alpha * p;

        LU_reverse(A, r, LAU);
        LAU = *_A.dot(LAU);
        LU_direct(A, LAU, LAU);

        beta = -(p * LAU) / pp;

        LU_reverse(A, r, U);

        z = U + beta * z;

        p = LAU + beta * p;
    }

    res.iters = k;
    res.residual = rr;

    r.clear();
    z.clear();
    p.clear();
    LAU.clear();
    U.clear();

    return res;
}

result solver::BCGSTAB_LU(matrix& _A, dvector& b, dvector& x)
{
    matrix A = _A;

    uint32_t dim = A.size();

    result res;

    dvector r0(dim), r(dim), z(dim), p(dim), y(dim);
    dvector LAUz(dim), LAUp(dim);

    double alpha, beta, gamma;

    x = _initial_approx;

    decomposer::LU(A);

    LU_direct(A, b - *_A.dot(x), r0);

    z = r0;
    r = r0;

    uint32_t k = 0;
    double rr = (r * r);

    for (; k < _max_iter && rr >= _eps; k++)
    {
        LU_reverse(A, z, LAUz);
        LAUz = *_A.dot(LAUz);
        LU_direct(A, LAUz, LAUz);

        double r_r0 = (r * r0);

        alpha = r_r0 / (LAUz * r0);

        p = r - alpha * LAUz;

        LU_reverse(A, p, LAUp);
        LAUp = *_A.dot(LAUp);
        LU_direct(A, LAUp, LAUp);

        gamma = (p * LAUp) / (LAUp * LAUp);
        y = y + alpha * z + gamma * p;

        r = p - gamma * LAUp;

        beta = alpha * (r * r0) / (gamma * r_r0);

        z = r + beta * z - beta * gamma * LAUz;

        LU_reverse(A, y, x);

        rr = (r * r);
    }

    res.iters = k;
    res.residual = rr;

    r.clear();
    r0.clear();
    z.clear();
    p.clear();
    y.clear();
    LAUp.clear();
    LAUz.clear();

    return res;
}

result solver::LOS(matrix& A, dvector& b, dvector& x)
{
    result res;

    x = _initial_approx;

    dvector r(A.size()), z(A.size()), p(A.size());

    r.resize(A.size());
    p.resize(A.size());
    z.resize(A.size());

    double alpha, beta;
    std::vector<double> Ar;
    r = b - (*A.dot(x));
    z = r;
    p = *A.dot(z);

    double rr = (r * r);

    uint32_t k = 0;

    for (; k < _max_iter && rr >= _eps; k++)
    {
        double pp = (p * p);
        alpha = (p * r) / pp;
        x = x + alpha * z;
        rr = (r * r) - alpha * alpha * pp;
        r = r - alpha * p;
        Ar = *A.dot(r);
        beta = -(p * Ar) / pp;
        z = r + beta * z;
        p = Ar + beta * p;
    }

    r.clear();
    p.clear();
    z.clear();

    res.iters = k;
    res.residual = rr;

    return res;
}