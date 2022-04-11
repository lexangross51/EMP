#include "decomposer.hpp"

void decomposer::LU(matrix& A)
{
    for (uint32_t i = 0; i < A.size(); i++)
    {
        double sum_di = 0.0;

        uint32_t i0 = A._ig[i];
        uint32_t i1 = A._ig[i + 1];

        for (uint32_t k = i0; k < i1; k++)
        {
            uint32_t j = A._jg[k];

            uint32_t j0 = A._ig[j];
            uint32_t j1 = A._ig[j + 1];

            uint32_t ik = i0;
            uint32_t kj = j0;

            double sum_l = 0.0;
            double sum_u = 0.0;

            while (ik < k && kj < j1)
            {
                if (A._jg[ik] == A._jg[kj])
                {
                    sum_l += A._ggl[ik] * A._ggu[kj];
                    sum_u += A._ggu[ik] * A._ggl[kj];
                    ik++;
                    kj++;
                }
                else
                    A._jg[ik] > A._jg[kj] ? kj++ : ik++;
            }

            A._ggl[k] -= sum_l;
            A._ggu[k] = (A._ggu[k] - sum_u) / A._di[j];
            sum_di += A._ggl[k] * A._ggu[k];
        }

        A._di[i] -= sum_di;
    }
}