#include "mfe.hpp"

mfe::mfe(space_grid* space, time_grid* time, function2D_t& u, function2D_t_r& f)
{
    _local_size = 9;

    _global_A = nullptr;
    _mesh_rz = nullptr;
    _mesh_rz = space;
    _time_layers = time;

    _local_G.resize(_local_size);
    _local_M.resize(_local_size);
    _local_b.resize(_local_size);

    for (uint8_t i = 0; i < _local_size; i++)
    {
        _local_G[i].resize(_local_size);
        _local_M[i].resize(_local_size);
    }

    _u = u;
    _f = f;

    rebuild_mesh(space);

    _q_prev1.resize(_mesh_rz->get_nodes_count());
    _q_prev2.resize(_mesh_rz->get_nodes_count());
    _q_prev3.resize(_mesh_rz->get_nodes_count());

    _global_A = new matrix(_mesh_rz->get_nodes_count());
    
    std::vector<uint32_t> ig, jg;

    portrait_generator::portrait(*_mesh_rz, ig, jg);

    _global_A->set_portrait(ig, jg);

    ig.clear();
    jg.clear();

    _global_b.resize(_mesh_rz->get_nodes_count());
}

void mfe::rebuild_mesh(space_grid* grid)
{
    _mesh_rz = grid;

    uint32_t nr = grid->get_width();
    uint32_t nz = grid->get_height();

    uint32_t n_nodes = (2 * nr - 1) * (2 * nz - 1);

    _mesh_rz = new space_grid((nr - 1) * (nz - 1), n_nodes);

    _mesh_rz->set_type(grid->get_type());
    _mesh_rz->set_width(2 * nr - 1);
    _mesh_rz->set_height(2 * nz - 1);

    // Формируем конечные элементы с новой нумерацией и сразу же 
    // добавляем координаты дополнительных точек
    std::vector<uint32_t> nodes(9);

    for (uint32_t ielem = 0; ielem < _mesh_rz->get_elems_count(); ielem++)
    {
        finite_element elem = grid->get_elem(ielem);

        // Формируем сам конечный элемент
        uint32_t global_num = 2 * uint32_t(ielem / (nr - 1)) * (2 * nr - 1) + 2 * (ielem % (nr - 1));

        for (uint8_t i = 0; i < 3; i++)
            nodes[i] = global_num + i;

        for (uint8_t i = 3, j = 0; i < 6; i++, j++)
            nodes[i] = global_num + 2 * nr - 1 + j;

        for (uint8_t i = 6, j = 0; i < 9; i++, j++)
            nodes[i] = global_num + 4 * nr - 2 + j;	

        _mesh_rz->add_elem(
            finite_element(nodes, elem.lambda, elem.sigma, elem.hi), ielem
        );

        // Формируем новые точки и заносим старые
        point2D p1 = grid->get_point(elem[0]);
        point2D p2 = grid->get_point(elem[1]);
        point2D p3 = grid->get_point(elem[2]);
        point2D p4 = grid->get_point(elem[3]);

        double r_aver = (p1.r + p2.r) / 2.0;
        double z_aver = (p1.z + p3.z) / 2.0;

        _mesh_rz->add_point(p1, nodes[0]);
        _mesh_rz->add_point(point2D(r_aver, p1.z), nodes[1]);
        _mesh_rz->add_point(p2, nodes[2]);
        _mesh_rz->add_point(point2D(p1.r, z_aver), nodes[3]);
        _mesh_rz->add_point(point2D(r_aver, z_aver), nodes[4]);
        _mesh_rz->add_point(point2D(p2.r, z_aver), nodes[5]);
        _mesh_rz->add_point(p3, nodes[6]);
        _mesh_rz->add_point(point2D(r_aver, p3.z), nodes[7]);
        _mesh_rz->add_point(p4, nodes[8]);
    }

    nodes.clear();


    // Меняем первые краевые
    std::vector<int> bc_nodes(grid->get_first_bound_count(), -1);
    std::vector<dirichlet_cond> dirichlet_nodes;

    for (uint32_t i = 0; i < grid->get_first_bound_count(); i++)
    {
        bool find = false;

        uint32_t bc_node = grid->get_dirichlet_cond(i).node;
        uint32_t new_bc_num;

        for (uint32_t ielem = 0; ielem < grid->get_elems_count() && !find; ielem++)
        {
            for (uint8_t node = 0; node < 4 && !find; node++)
            {
                if (grid->get_elem(ielem)[node] == bc_node)
                {
                    if (node == 0) new_bc_num = _mesh_rz->get_elem(ielem)[0];
                    if (node == 1) new_bc_num = _mesh_rz->get_elem(ielem)[2];
                    if (node == 2) new_bc_num = _mesh_rz->get_elem(ielem)[6];
                    if (node == 3) new_bc_num = _mesh_rz->get_elem(ielem)[8];

                    bc_nodes[i] = new_bc_num;

                    find = true;
                }
            }
        }
    }

    if (bc_nodes.size())
    {
        dirichlet_nodes.resize(2 * bc_nodes.size());

        uint32_t bc_num = 0;

        for (const auto& node : bc_nodes)
            dirichlet_nodes[2 * bc_num++] = { (uint32_t)node, 0.0 };

        for (uint32_t i = 0; i < bc_nodes.size() - 1; i++)
            dirichlet_nodes[2 * i + 1] = { 
                (dirichlet_nodes[2 * i].node + dirichlet_nodes[2 * i + 2].node) / 2, 0.0 
            };

        if ( 2 * (nr + nz) - 4 == grid->get_first_bound_count() && grid->get_second_bound_count() == 0)
            dirichlet_nodes[dirichlet_nodes.size() - 1] = {
                (dirichlet_nodes[0].node + dirichlet_nodes[dirichlet_nodes.size() - 2].node) / 2, 0.0
            };

        _mesh_rz->set_dirichlet_conds(dirichlet_nodes);

        dirichlet_nodes.clear();
        bc_nodes.clear();
    }

    // Копируем вторые краевые
    if (grid->get_second_bound_count())
    {
        std::vector<neumann_cond> neumann_nodes(grid->get_second_bound_count());

        for (uint32_t i = 0; i < neumann_nodes.size(); i++)
        {
            neumann_nodes[i] = grid->get_neumann_cond(i);
            neumann_nodes[i].local_node_2 = 6;
        }

        _mesh_rz->set_neumann_conds(neumann_nodes);

        neumann_nodes.clear();
    }

    //grid->~space_grid();
}

void mfe::build_local_matrices(finite_element& elem)
{
    omega omega = { _mesh_rz->get_point(elem[0]), _mesh_rz->get_point(elem[_local_size - 1]) };

    basis_function bf(omega);

    function2D f;

    for (uint8_t i = 0; i < _local_size; i++)
    {
        for (uint8_t j = 0; j < _local_size; j++)
        {
            f = [&](double r, double z)
            {
                point2D point(r, z);

                double dpsi_i_r = bf.d_psi(i, 1, point);
                double dpsi_j_r = bf.d_psi(j, 1, point);

                double dpsi_i_z = bf.d_psi(i, 2, point);
                double dpsi_j_z = bf.d_psi(j, 2, point);

                return (dpsi_i_r * dpsi_j_r + dpsi_i_z * dpsi_j_z) * r;
            };

            _local_G[i][j] = _local_G[j][i] = _gauss.integrate2D(f, omega);
        }
    }

    for (uint8_t i = 0; i < _local_size; i++)
    {
        for (uint8_t j = 0; j < _local_size; j++)
        {
            f = [&](double r, double z)
            {
                point2D point(r, z);

                double psi_i = bf.psi(i, point);
                double psi_j = bf.psi(j, point);

                return psi_i * psi_j * r;
            };

            _local_M[i][j] = _local_M[j][i] = _gauss.integrate2D(f, omega);
        }
    }
}

void mfe::build_local_vector(finite_element& elem, double time)
{
    dvector f(_local_size);

    for (uint32_t i = 0; i < _local_b.size(); i++)
        _local_b[i] = 0.0;

    for (uint8_t i = 0; i < _local_size; i++)
    {
        point2D point = _mesh_rz->get_point(elem[i]);

        f[i] = _f(point.r, point.z, time, elem.lambda, elem.sigma, elem.hi);
    }

    for (uint32_t i = 0; i < _local_size; i++)
        for (uint32_t j = 0; j < _local_size; j++)
            _local_b[i] += _local_M[i][j] * f[j];

    f.clear();
}

void mfe::add_to_global_matrix(const uint32_t i, const uint32_t j, const double val)
{
    if (i == j)
    {
        _global_A->di_a(i, val);
        return;
    }

    if (i < j)
    {
        for (uint32_t ind = _global_A->ig(j); ind < _global_A->ig(j + 1); ind++)
        {
            if (_global_A->jg(ind) == i)
            {
                _global_A->ggu_a(ind, val);
                return;
            }
        }
    }
    else
    {
        for (uint32_t ind = _global_A->ig(i); ind < _global_A->ig(i + 1); ind++)
        {
            if (_global_A->jg(ind) == j)
            {
                _global_A->ggl_a(ind, val);
                return;
            }
        }
    }
}

void mfe::initial_condition()
{
    double t0 = _time_layers->get_time(0);
    double t1 = _time_layers->get_time(1);
    double t2 = _time_layers->get_time(2);

    for (uint32_t i = 0; i < _mesh_rz->get_nodes_count(); i++)
    {
        point2D point = _mesh_rz->get_point(i);

        _q_prev1[i] = _u(point.r, point.z, t2);
        _q_prev2[i] = _u(point.r, point.z, t1);
        _q_prev3[i] = _u(point.r, point.z, t0);
    }
}

void mfe::assembly_global_matrix_and_vector(const uint32_t time_m)
{
    // Зануляем матрицу и вектор
    _global_A->null();

    for (uint32_t i = 0; i < _global_b.size(); i++)
        _global_b[i] = 0.0;

    // Собираем глобальную матрицу и сразу делаем добавки от четырехслойной схемы
    double t = _time_layers->get_time(time_m);
    double t1 = _time_layers->get_time(time_m - 1);
    double t2 = _time_layers->get_time(time_m - 2);
    double t3 = _time_layers->get_time(time_m - 3);

    double d_t03 = t - t3;
    double d_t02 = t - t2;
    double d_t01 = t - t1;
    double d_t13 = t1 - t3;
    double d_t12 = t1 - t2;
    double d_t23 = t2 - t3;

    double M_q_prev3 = 0;
    double M_q_prev2 = 0;
    double M_q_prev1 = 0;

    for (uint32_t ielem = 0; ielem < _mesh_rz->get_elems_count(); ielem++)
    {
        finite_element elem = _mesh_rz->get_elem(ielem);

        build_local_matrices(elem);

        build_local_vector(elem, t);

        double additive_M = elem.sigma * (1.0 / d_t03 + 1.0 / d_t02 + 1.0 / d_t01) +
                        2 * elem.hi * (d_t01 + d_t02 + d_t03) / (d_t01 * d_t02 * d_t03);

        M_q_prev3 = (elem.sigma * (d_t02 * d_t01) + 2 * elem.hi * (d_t02 + d_t01)) / (d_t23 * d_t13 * d_t03);
        M_q_prev2 = (elem.sigma * (d_t03 * d_t01) + 2 * elem.hi * (d_t03 + d_t01)) / (d_t23 * d_t12 * d_t02);
        M_q_prev1 = (elem.sigma * (d_t03 * d_t02) + 2 * elem.hi * (d_t03 + d_t02)) / (d_t13 * d_t12 * d_t01);

        for (uint32_t i = 0; i < _local_size; i++)
        {
            _global_b[elem[i]] += _local_b[i];

            for (uint32_t j = 0; j < _local_size; j++)
            {
                add_to_global_matrix(
                    elem[i], elem[j], 
                    elem.lambda * _local_G[i][j] + additive_M * _local_M[i][j]
                );

                _global_b[elem[i]] += M_q_prev3 * _local_M[i][j] * _q_prev3[elem[j]]
                                    - M_q_prev2 * _local_M[i][j] * _q_prev2[elem[j]]
                                    + M_q_prev1 * _local_M[i][j] * _q_prev1[elem[j]];
            }
        }
    }
}

void mfe::add_dirichlet(const uint32_t time_m)
{
    for (uint32_t i = 0; i < _mesh_rz->get_first_bound_count(); i++)
    {
        uint32_t bc_node = _mesh_rz->get_dirichlet_cond(i).node;

        point2D point = _mesh_rz->get_point(bc_node);

        // На диагональ 1, в правую часть значение функции
        _global_A->di_s(bc_node, 1.0);
        _global_b[bc_node] = _u(point.r, point.z, _time_layers->get_time(time_m));

        // Обнуляем строку
        uint32_t i0 = _global_A->ig(bc_node);
        uint32_t i1 = _global_A->ig(bc_node + 1);

        for (uint32_t k = i0; k < i1; k++)
            _global_A->ggl_s(k, 0.0);

        for (uint32_t k = bc_node + 1; k < _global_A->size(); k++)
        {
            uint32_t i0 = _global_A->ig(k);
            uint32_t i1 = _global_A->ig(k + 1);

            for (uint32_t j = i0; j < i1; j++)
            {
                if (_global_A->jg(j) == bc_node)
                    _global_A->ggu_s(j, 0.0);
            }
        }
    }
}

void mfe::add_meumann(const uint32_t time_m)
{
    function2D f;

    double time = _time_layers->get_time(time_m);

    if (_mesh_rz->get_second_bound_count())
    {
        for (uint32_t i = 0; i < _mesh_rz->get_second_bound_count(); i++)
        {
            auto bcond = _mesh_rz->get_neumann_cond(i);

            auto elem = _mesh_rz->get_elem(bcond.ielem);
            uint32_t loc_node1 = bcond.local_node_1;
            uint32_t loc_node2 = bcond.local_node_2;
            double tetta = bcond.tetta;

            uint32_t ind_1 = elem[loc_node1];
            uint32_t ind_2 = elem[loc_node2];

            double rk = _mesh_rz->get_point(elem[0]).r;
            double rk1 = _mesh_rz->get_point(elem[8]).r;
            double hr = rk1 - rk;

            double zk = _mesh_rz->get_point(elem[0]).z;
            double zk1 = _mesh_rz->get_point(elem[6]).z;
            double z_aver = (zk + zk1) / 2.0;

            omega omega = { point2D(rk, zk), point2D(rk, zk1) };


            double du_dr = (_u(rk + hr, zk, time) - _u(rk - hr, zk, time)) / (2 * hr);
            tetta *= du_dr;

            f = [=](double r, double z)
            {
                return 2.0 * (z - z_aver) * (z - zk1) / ((zk1 - zk) * (zk1 - zk));
            };
            _global_b[ind_1] += tetta * rk * _gauss.integrate1D(f, omega);

            f = [=](double r, double z)
            {
                return 2.0 * (z - zk) * (z - z_aver) / ((zk1 - zk) * (zk1 - zk));
            };
            _global_b[ind_2] += tetta * rk * _gauss.integrate1D(f, omega);

            f = [=](double r, double z)
            {
                return -4.0 * (z - zk) * (z - zk1) / ((zk1 - zk) * (zk1 - zk));
            };
            _global_b[(ind_1 + ind_2) / 2] += tetta * rk * _gauss.integrate1D(f, omega);
        }
    }
}

double mfe::solve()
{
    solver slv;

    initial_condition();

    for (uint32_t t = 3; t < _time_layers->get_nodes_count(); t++)
    {
        calc_exact(t);

        assembly_global_matrix_and_vector(t);

        add_meumann(t);

        add_dirichlet(t);

        slv.init(100, 1e-17, _q_prev1, iterative_method::LOS_LU);
        auto res = slv.solve(*_global_A, _global_b, _q);

        _q_prev3 = _q_prev2;
        _q_prev2 = _q_prev1;
        _q_prev1 = _q;
    }

    return error();
}

void mfe::calc_exact(const uint32_t time_m)
{
    _exact.resize(_mesh_rz->get_nodes_count());

    double t = _time_layers->get_time(time_m);

    for (uint32_t i = 0; i < _exact.size(); i++)
    {
        point2D point = _mesh_rz->get_point(i);

        _exact[i] = _u(point.r, point.z, t);
    }
}

double mfe::error()
{
    auto diff = _exact - _q;

    return norm(diff) / norm(_exact);
}