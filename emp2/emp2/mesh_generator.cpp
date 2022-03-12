#include "mesh_generator.h"

void mesh_generator::input(std::string dir)
{
    std::ifstream in(dir + "data.txt");

    if (in.is_open())
    {
        in >> area_start >> area_end >> sigma;
        in >> nx >> kx;
        in >> area_nested;

        in >> time_start >> time_end;
        in >> nt >> kt;
        in >> time_nested;

        in.close();
    }
    else
        throw "Can't open file";
}

void mesh_generator::generate_XT(mesh::mesh_type area_type, mesh::mesh_type time_type)
{
    // Делаем разбиения по пространству
    if (area_nested == 1)
    {
        nx *= 2;
        kx = sqrt(kx);
    }
    else if (area_nested == 2)
    {
        nx *= 4;
        kx = sqrt(sqrt(kx));
    }
    else if (area_nested == 3)
    {
        nx *= 8;
        kx = sqrt(sqrt(sqrt(kx)));
    }

    X.resize(nx + 1);
    if (area_type == mesh::mesh_type::UNIFORM)
    {
        double hx = (area_end - area_start) / nx;

        X[0] = area_start;
        for (uint32_t i = 1; i < nx + 1; i++)
        {
            area_start += hx;
            X[i] = area_start;
        }
    }
    else
    {
        double hx = (area_end - area_start) * (1 - kx) / (1 - pow(kx, nx));

        X[0] = area_start;

        for (uint32_t i = 1; i < nx + 1; i++)
        {
            area_start += hx;
            X[i] = area_start;
            hx *= kx;
        }
    }

    // Делаем разбиения по времени
    if (time_nested == 1)
    {
        nt *= 2;
        kt = sqrt(kt);
    }
    else if (time_nested == 2)
    {
        nt *= 4;
        kt = sqrt(sqrt(kt));
    }
    else if (time_nested == 3)
    {
        nt *= 8;
        kt = sqrt(sqrt(sqrt(kt)));
    }

    T.resize(nt + 1);
    if (time_type == mesh::mesh_type::UNIFORM)
    {
        double ht = (time_end - time_start) / nt;

        T[0] = time_start;
        for (uint32_t i = 1; i < nt + 1; i++)
        {
            time_start += ht;
            T[i] = time_start;
        }
    }
    else
    {
        double ht = (time_end - time_start) * (1 - kt) / (1 - pow(kt, nt));

        T[0] = time_start;

        for (uint32_t i = 1; i < nt + 1; i++)
        {
            time_start += ht;
            T[i] = time_start;
            ht *= kt;
        }
    }
}

void mesh_generator::build_mesh(mesh& mesh, mesh::mesh_type area_type, mesh::mesh_type time_type)
{
    input();

    generate_XT(area_type, time_type);

    mesh.set_time(T);

    for (uint32_t i = 0; i < X.size() - 1; i++)
        mesh.add_elem(finite_elem(i, X[i], X[i + 1], sigma));

    mesh.save();
}