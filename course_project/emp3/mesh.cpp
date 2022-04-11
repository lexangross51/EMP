#include "mesh.hpp"

// =========================================================
void space_grid::save(std::string path)
{
    std::ofstream grid_out(path + "mesh.json");
    nlohmann::json grid{};

    if (grid_out.is_open())
    {
        for (uint32_t i = 0; i < _elems.size(); i++)
            grid["elems"].push_back(_elems[i].nodes);

        for (uint32_t i = 0; i < _points.size(); i++)
            grid["points"].push_back({ _points[i].r, _points[i].z });

        grid_out << grid << std::endl;

        grid_out.close();
    }
    else throw "Can't open file\n";
}

// =========================================================
void time_grid::save(std::string path)
{
    std::ofstream time_out(path + "time.json");
    nlohmann::json time{};

    if (time_out.is_open())
    {
        time["time"] = _layers;

        time_out << time;
        
        time_out.close();
    }
    else throw "Can't open file\n";
}