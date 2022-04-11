#include "tests.hpp"

tests::tests()
{
    _u_rz.resize(6);
    _u_rz_names.resize(6);
    _u_t.resize(6);
    _u_t_names.resize(6);
    _f.resize(6);

    for (auto& it : _f) it.resize(6);

    _u_rz_names[0] = "1";
    _u_rz_names[1] = "z";
    _u_rz_names[2] = "r^2 - 2z^2";
    _u_rz_names[3] = "r^3 + z^3";
    _u_rz_names[4] = "r^4 + z^3";
    _u_rz_names[5] = "e^(r^2) + z";


    _u_rz[0] = [](double r, double z) { return 1; };
    _u_rz[1] = [](double r, double z) { return z; };
    _u_rz[2] = [](double r, double z) { return r*r - 2*z*z; };
    _u_rz[3] = [](double r, double z) { return 2*r*r + 3*z*z*z; };
    _u_rz[4] = [](double r, double z) { return r*r*r*r + z*z*z; };
    _u_rz[5] = [](double r, double z) { return exp(r*r) + z; };


    _u_t_names[0] = "1";
    _u_t_names[1] = "t";
    _u_t_names[2] = "t^2";
    _u_t_names[3] = "t^3";
    _u_t_names[4] = "t^4";
    _u_t_names[5] = "sin(t)";

    _u_t[0] = [](double t) { return 1; };
    _u_t[1] = [](double t) { return t; };
    _u_t[2] = [](double t) { return t*t; };
    _u_t[3] = [](double t) { return t*t*t; };
    _u_t[4] = [](double t) { return t*t*t*t; };
    _u_t[5] = [](double t) { return sin(t); };


    _f[0][0] = [](double r, double z, double t, double lambda, double sigma, double hi) {return 0; };
    _f[0][1] = [](double r, double z, double t, double lambda, double sigma, double hi) {return 0; };
    _f[0][2] = [](double r, double z, double t, double lambda, double sigma, double hi) {return 0; };
    _f[0][3] = [](double r, double z, double t, double lambda, double sigma, double hi) {return -lambda * (8 + 18*z); };
    _f[0][4] = [](double r, double z, double t, double lambda, double sigma, double hi) {return -lambda * (16*r*r + 6*z); };
    _f[0][5] = [](double r, double z, double t, double lambda, double sigma, double hi) {return -lambda * 4*exp(r*r) * (1 + r*r); };

    _f[1][0] = [](double r, double z, double t, double lambda, double sigma, double hi) {return sigma; };
    _f[1][1] = [](double r, double z, double t, double lambda, double sigma, double hi) {return sigma; };
    _f[1][2] = [](double r, double z, double t, double lambda, double sigma, double hi) {return sigma; };
    _f[1][3] = [](double r, double z, double t, double lambda, double sigma, double hi) {return sigma - lambda * (8 + 18 * z); };
    _f[1][4] = [](double r, double z, double t, double lambda, double sigma, double hi) {return sigma - lambda * (16 * r * r + 6 * z); };
    _f[1][5] = [](double r, double z, double t, double lambda, double sigma, double hi) {return sigma - lambda * 4 * exp(r * r) * (1 + r * r); };

    _f[2][0] = [](double r, double z, double t, double lambda, double sigma, double hi) {return 2*t*sigma + 2*hi; };
    _f[2][1] = [](double r, double z, double t, double lambda, double sigma, double hi) {return 2*t*sigma + 2*hi; };
    _f[2][2] = [](double r, double z, double t, double lambda, double sigma, double hi) {return 2*t*sigma + 2*hi; };
    _f[2][3] = [](double r, double z, double t, double lambda, double sigma, double hi) {return 2*t*sigma + 2*hi - lambda * (8 + 18 * z); };
    _f[2][4] = [](double r, double z, double t, double lambda, double sigma, double hi) {return 2*t*sigma + 2*hi - lambda * (16 * r * r + 6 * z); };
    _f[2][5] = [](double r, double z, double t, double lambda, double sigma, double hi) {return 2*t*sigma + 2*hi - lambda * 4 * exp(r * r) * (1 + r * r); };

    _f[3][0] = [](double r, double z, double t, double lambda, double sigma, double hi) {return 3*t*t * sigma + 6*t * hi; };
    _f[3][1] = [](double r, double z, double t, double lambda, double sigma, double hi) {return 3*t*t * sigma + 6*t * hi; };
    _f[3][2] = [](double r, double z, double t, double lambda, double sigma, double hi) {return 3*t*t * sigma + 6*t * hi; };
    _f[3][3] = [](double r, double z, double t, double lambda, double sigma, double hi) {return 3*t*t * sigma + 6*t * hi - lambda * (8 + 18 * z); };
    _f[3][4] = [](double r, double z, double t, double lambda, double sigma, double hi) {return 3*t*t * sigma + 6*t * hi - lambda * (16 * r * r + 6 * z); };
    _f[3][5] = [](double r, double z, double t, double lambda, double sigma, double hi) {return 3*t*t * sigma + 6*t * hi - lambda * 4 * exp(r * r) * (1 + r * r); };

    _f[4][0] = [](double r, double z, double t, double lambda, double sigma, double hi) {return 4*t*t*t * sigma + 12*t*t * hi; };
    _f[4][1] = [](double r, double z, double t, double lambda, double sigma, double hi) {return 4*t*t*t * sigma + 12*t*t * hi; };
    _f[4][2] = [](double r, double z, double t, double lambda, double sigma, double hi) {return 4*t*t*t * sigma + 12*t*t * hi; };
    _f[4][3] = [](double r, double z, double t, double lambda, double sigma, double hi) {return 4*t*t*t * sigma + 12*t*t * hi - lambda * (8 + 18 * z); };
    _f[4][4] = [](double r, double z, double t, double lambda, double sigma, double hi) {return 4*t*t*t * sigma + 12*t*t * hi - lambda * (16 * r * r + 6 * z); };
    _f[4][5] = [](double r, double z, double t, double lambda, double sigma, double hi) {return 4*t*t*t * sigma + 12*t*t * hi - lambda * 4 * exp(r * r) * (1 + r * r); };

    _f[5][0] = [](double r, double z, double t, double lambda, double sigma, double hi) {return cos(t) * sigma - sin(t) * hi; };
    _f[5][1] = [](double r, double z, double t, double lambda, double sigma, double hi) {return cos(t) * sigma - sin(t) * hi; };
    _f[5][2] = [](double r, double z, double t, double lambda, double sigma, double hi) {return cos(t) * sigma - sin(t) * hi; };
    _f[5][3] = [](double r, double z, double t, double lambda, double sigma, double hi) {return cos(t) * sigma - sin(t) * hi - lambda * (8 + 18 * z); };
    _f[5][4] = [](double r, double z, double t, double lambda, double sigma, double hi) {return cos(t) * sigma - sin(t) * hi - lambda * (16 * r * r + 6 * z); };
    _f[5][5] = [](double r, double z, double t, double lambda, double sigma, double hi) {return cos(t) * sigma - sin(t) * hi - lambda * 4 * exp(r * r) * (1 + r * r); };
}

void tests::run()
{
    std::ofstream file(directory + "tests\\table.txt", std::ios::app);

    time_grid_generator tgg;
    time_grid* tg = nullptr;

    space_grid_generator sgg;
    space_grid* sg = nullptr;

    tgg.build_mesh(tg);

    sgg.build_mesh(sg);
    sg->save();

    file << std::left << std::setw(12) << "time/space";
    for (uint8_t i = 0; i < _u_rz_names.size(); i++)
        file << "|" << std::setw(12) << _u_rz_names[i];
    file << std::endl;
    file << "---------------------------------------------";
    file << "---------------------------------------------";
    file << std::endl;


    for (uint8_t i = 0; i < _u_t.size(); i++)
    {
        file << std::left << std::setw(12) << _u_t_names[i];

        for (uint8_t j = 0; j < _u_rz.size(); j++)
        {
            function2D_t u = [&](double r, double z, double t) { 
                return _u_rz[j](r, z) + _u_t[i](t); 
            };

            mfe mfe(sg, tg, u, _f[i][j]);
            double error = mfe.solve();

            file << "|" << std::left << std::setw(12) << error;
        }

        file << std::endl;
    }

    file << std::endl << std::endl << std::endl;
    file.close();
}