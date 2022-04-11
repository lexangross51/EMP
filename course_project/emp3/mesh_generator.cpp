#include "mesh_generator.hpp"

// ===================================================================
space_grid_generator::space_grid_generator()
{
    _rw = _rb = 0.0;
    _depth = 0.0;
    _nr = 0;
    _kr = 0.0;
    _type = mesh_type::UNIFORM;
    _nested = 0;
}

void space_grid_generator::read_data(std::string path)
{
    std::ifstream input(path + "space_grid.json");

    if (input.is_open())
    {
        nlohmann::json space_grid{};

        input >> space_grid;

        input.close();

        _rw		= space_grid["area"]["Rw"];
        _rb		= space_grid["area"]["Rb"];
        _depth	= space_grid["area"]["depth"];

        auto height		= space_grid["area"]["heights"];
        auto lambda		= space_grid["area"]["lambda"];
        auto sigma		= space_grid["area"]["sigma"];
        auto hi			= space_grid["area"]["hi"];

        _nr				= space_grid["parameters"]["nr"];
        _kr				= space_grid["parameters"]["kr"];
        auto z_parts	= space_grid["parameters"]["nz"];
        auto z_coeffs	= space_grid["parameters"]["kz"];

        _heights.resize(height.size());
        _lambda_s.resize(height.size());
        _sigma_s.resize(height.size());
        _hi_s.resize(height.size());
        _nz.resize(height.size());
        _kz.resize(height.size());

        for (uint16_t i = 0; i < _heights.size(); i++)
        {
            _heights[i]	 = height[i];
            _lambda_s[i] = lambda[i];
            _sigma_s[i]	 = sigma[i];
            _hi_s[i]	 = hi[i];
            _nz[i]		 = z_parts[i];
            _kz[i]		 = z_coeffs[i];
        }

        _type	= space_grid["parameters"]["type"];
        _nested	= space_grid["parameters"]["nested"];
        
        // Краевые
        _bconds[0] = { 1, space_grid["boundary"]["border_1"], {} };
        _bconds[1] = { 2, space_grid["boundary"]["border_2"], {} };
        _bconds[2] = { 3, space_grid["boundary"]["border_3"], {} };

        auto bcond = space_grid["boundary"]["border_4"];

        if (bcond.size() > 1)
        {
            boundary_cond bc;
            bc.border = 4;
            bc.type = boundary_type(2);

            bc.boundaries.resize(bcond.size());

            for (uint32_t i = 0; i < bcond.size() / 2; i++)
            {
                bc.boundaries[2 * i] = bcond[2 * i];
                bc.boundaries[2 * i + 1] = bcond[2 * i + 1];
            }

            _bconds[3] = bc;
        }
        else
        {
            boundary_cond bc;
            bc.border = 4;
            bc.type = boundary_type(1);

            _bconds[3] = bc;
        }
    }
    else throw "Can't open file space_grid.json\n";
}

void space_grid_generator::generate_nodes()
{
    // Если вложенная, то увеличиваем число разбиений
    if (_nested == 1)
    {
        _nr *= 2;
        _kr = sqrt(_kr);
        for (auto& it : _nz) it *= 2;
        for (auto& it : _kz) it = sqrt(it);
    }
    else if (_nested == 2)
    {
        _nr *= 4;
        _kr = sqrt(sqrt(_kr));
        for (auto& it : _nz) it *= 4;
        for (auto& it : _kz) it = sqrt(sqrt(it));
    }
    else if (_nested == 3)
    {
        _nr *= 8;
        _kr = sqrt(sqrt(sqrt(_kr)));
        for (auto& it : _nz) it *= 8;
        for (auto& it : _kz) it = sqrt(sqrt(sqrt(it)));
    }

    // Если сетка равномерная
    if (_type == mesh_type::UNIFORM)
    {
        double height = 0.0;

        double z_start = _depth;

        for (const auto& it : _heights) height += it;

        double hr = (_rb - _rw) / _nr;
        double hz = height / _nz[0];

        for (uint32_t i = 0; i < _nr + 1; i++)
            _r.insert(_rw + hr * i);

        for (uint32_t i = 0; i < _nz[0] + 1; i++)
            _z.insert(z_start + i * hz);
    }
    // Если сетка неравномерная
    else
    {
        double hr = 0.0, hz = 0.0;

        if (_kr == 1)
            hr = (_rb - _rw) / _nr;
        else
            hr = (_rb - _rw) * (1 - _kr) / (1 - pow(_kr, _nr));
        
        double r_start = _rw;

        for (uint32_t i = 0; i < _nr + 1; i++)
        {
            _r.insert(r_start);
            r_start += hr;
            hr *= _kr;
        }


        double z_start = _depth;

        for (uint32_t i = 0; i < _heights.size(); i++)
        {
            if (_kz[i] == 1)
                hz = _heights[i] / _nz[i];
            else
                hz = (_heights[i]) * (1 - _kz[i]) / (1 - pow(_kz[i], _nz[i]));

            double coord = z_start;

            for (uint32_t j = 0; j < _nz[i] + 1; j++)
            {
                _z.insert(coord);
                coord += hz;
                hz *= _kz[i];
            }

            z_start += _heights[i];
        }
    }

    // Добавляем положениия зоны перфорации
    for (uint8_t i = 0; i < 4; i++)
    {
        auto bcond = _bconds[i];

        if (bcond.type == boundary_type::NEUMANN)
        {
            for (uint32_t j = 0; j < bcond.boundaries.size() / 2; j++)
            {
                _z.insert(bcond.boundaries[2 * j]);
                _z.insert(bcond.boundaries[2 * j + 1]);
            }
        }
    }
}

void space_grid_generator::make_bc_nodes(space_grid& grid)
{
    std::vector<uint32_t> nodes;
    std::vector<dirichlet_cond> dirichlet_nodes;

    std::vector<neumann_cond> neumann_nodes;
    std::vector<uint32_t> elems;

    for (const auto& bc : _bconds)
    {
        if (bc.type == boundary_type::DIRICHLET)
        {
            switch (bc.border)
            {
            case 1: {
                for (uint32_t i = 0; i < _nr + 1; i++)
                    nodes.push_back(i);
                break;
            }

            case 2: {
                for (uint32_t i = 0; i < _z.size(); i++)
                    nodes.push_back(_nr + i * (_nr + 1));
                break;
            }

            case 3: {
                for (int i = _nr; i >= 0; i--)
                    nodes.push_back( (_z.size() - 1) * (_nr + 1) + i);
                break;
            }
                  
            case 4: {
                for (int i = _z.size() - 1; i > 0 ; i--)
                    nodes.push_back(i * (_nr + 1));
                break;
            }

            default:
                throw "Undefined border!!!\n";
            }
        }
        else if (bc.type == boundary_type::NEUMANN)
        {
            for (uint32_t i = 0; i < bc.boundaries.size() / 2; i++)
            {
                uint32_t elem_start = 0;
                uint32_t elem_end = 0;

                uint32_t finded_pos = 0;

                for (const auto& it : _z)
                {
                    if (abs(bc.boundaries[2 * i] - it) < 1e-14)
                        elem_start = finded_pos;
                    else if (abs(bc.boundaries[2 * i + 1] - it) < 1e-14)
                        elem_end = finded_pos;

                    finded_pos++;
                }

                elems.push_back(elem_start++ * _nr);

                for (; elem_start < elem_end; elem_start++)
                    elems.push_back(elem_start++ * _nr);
            }
        }
    }

    if (nodes.size())
    {
        nodes.erase(std::unique(nodes.begin(), nodes.end()), nodes.end());

        dirichlet_nodes.resize(nodes.size());

        uint32_t bc1_num = 0;

        for (const auto& node : nodes)
            dirichlet_nodes[bc1_num++] = { node, 0.0 };

        grid.set_dirichlet_conds(dirichlet_nodes);

        nodes.clear();
        dirichlet_nodes.clear();
    }

    if (elems.size())
    {
        neumann_nodes.resize(elems.size());

        for (uint32_t i = 0; i < neumann_nodes.size(); i++)
            neumann_nodes[i] = { elems[i], 0, 2, -1.0 };

        grid.set_neumann_conds(neumann_nodes);

        elems.clear();
        neumann_nodes.clear();
    }
}

void space_grid_generator::build_mesh(space_grid*& grid)
{
    read_data();
    generate_nodes();

    grid = new space_grid(_nr * (_z.size() - 1), _r.size() * _z.size());
    
    grid->set_width(_nr + 1);
    grid->set_height(_z.size());
    grid->set_type(_type);

    // Формируем узлы
    uint32_t node = 0;

    for (const auto& z_p : _z)
        for (const auto& r_p : _r)
            grid->add_point(point2D(r_p, z_p), node++);

    // Формируем конечные элементы
    uint32_t ielem = 0;

    std::vector<uint32_t> nodes(4);

    uint32_t i = 0;

    auto z = _z.begin();

    for (uint32_t layer = 0; layer < _heights.size(); layer++)
    {
        for (; abs(z._Ptr->_Myval - (_depth + _heights[layer])) >= 1e-14; z++, i++)
        {
            for (uint32_t j = 0; j < _nr; j++)
            {
                nodes[0] = j + (_nr + 1) * i;
                nodes[1] = j + (_nr + 1) * i + 1;
                nodes[2] = j + (_nr + 1) * i + _nr + 1;
                nodes[3] = j + (_nr + 1) * i + _nr + 2;

                grid->add_elem(finite_element(
                    nodes, _lambda_s[layer], _sigma_s[layer], _hi_s[layer]
                ), ielem++);
            }
        }
        _depth += _heights[layer];
    }

    _lambda_s.clear();
    _sigma_s.clear();
    _hi_s.clear();

    _kz.clear();
    _nz.clear();

    make_bc_nodes(*grid);

    _r.clear();
    _z.clear();
    _heights.clear();
}

// ===================================================================
time_grid_generator::time_grid_generator()
{
    _time_begin = _time_end = 0.0;
    _nt = 0;
    _kt = 0.0;
    _type = mesh_type::UNIFORM;
    _nested = 0;
}

void time_grid_generator::read_data(std::string path)
{
    std::ifstream input(path + "time_grid.json");

    if (input.is_open())
    {
        nlohmann::json time{};

        input >> time;

        _time_begin	= time["time_interval"]["time_begin"];
        _time_end	= time["time_interval"]["time_end"];
        
        _nt			= time["parameters"]["nt"];
        _kt			= time["parameters"]["kt"];
        _type		= time["parameters"]["type"];
        _nested		= time["parameters"]["nested"];

        input.close();
    }
    else throw "Can't opne file time_grid.json\n";
}

void time_grid_generator::generate_nodes()
{
    if (_nested == 1) {
        _nt *= 2; _kt = sqrt(_kt);
    }
    else if (_nested == 2) {
        _nt *= 4; _kt = sqrt(sqrt(_kt));
    }
    else if (_nested == 3) {
        _nt *= 8; _kt = sqrt(sqrt(sqrt(_kt)));
    }

    if (_type == mesh_type::UNIFORM)
    {
        double ht = (_time_end - _time_begin) / _nt;

        _time_layers.resize(_nt + 1);

        for (uint32_t i = 0; i < _time_layers.size(); i++)
            _time_layers[i] = _time_begin + ht * i;
    }
    else
    {
        double ht = (_time_end - _time_begin) * (1 - _kt) / (1 - pow(_kt, _nt));

        _time_layers.resize(_nt + 1);

        for (uint32_t i = 0; i < _time_layers.size(); i++)
            _time_layers[i] = _time_begin + ht * _kt * i;
    }
}

void time_grid_generator::build_mesh(time_grid*& grid)
{
    read_data();
    generate_nodes();

    grid = new time_grid(_nt + 1);
    grid->_layers = _time_layers;
    grid->set_type(_type);

    _time_layers.clear();
}