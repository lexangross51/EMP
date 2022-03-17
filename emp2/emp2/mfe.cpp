#include "mfe.h"

mfe::mfe(mesh& mesh)
{
    global_A = new matrix(mesh.get_fe_count() + 1, 3);
    A = nullptr;
    global_f.resize(mesh.get_fe_count() + 1);
    q.resize(mesh.get_fe_count() + 1);
    q_prev.resize(mesh.get_fe_count() + 1);

    method_to_solve = method::SIMPLE_ITERATION;

    local_f.resize(2);
    local_A.resize(2);
    for (auto& row : local_A)
        row.resize(2);
}

void mfe::set_functions(const function& _u, const function_f& _f, const function_lambda& _lambda)
{
    u = _u;
    f = _f;
    lambda = _lambda;
}

// Сборка локальной матрицы
void mfe::build_local_matrix(const finite_elem& elem, const double delta_t)
{
    double h = elem.x_next - elem.x;

    //double du_dx = (q_prev[elem.number + 1] - q_prev[elem.number]) / h;


    double lambda1 = lambda(q_prev[elem.number + 1] / h, q_prev[elem.number] / h, elem.x);
    double lambda2 = lambda(q_prev[elem.number + 1] / h, q_prev[elem.number] / h, elem.x_next);

    // Матрица жесткости
    double coef = (lambda1 + lambda2) / (2 * h);
    local_A[0][0] = local_A[1][1] = coef;
    local_A[0][1] = local_A[1][0] = -coef;

    // Матрица жесткости + матрица масс
    coef = elem.sigma * h / (6 * delta_t);
    local_A[0][0] += 2 * coef;
    local_A[1][1] += 2 * coef;
    local_A[0][1] += coef;
    local_A[1][0] += coef;
}

// Сборка локального вектора
void mfe::build_local_vector(const finite_elem& elem, const double t, const double delta_t)
{
    double h = elem.x_next - elem.x;

    double f1 = f(elem.x, t, elem.sigma);
    double f2 = f(elem.x_next, t, elem.sigma);

    local_f[0] = h * (2 * f1 + f2) / 6.0 + elem.sigma * h * (2 * q_prev[elem.number] + q_prev[elem.number + 1]) / (6 * delta_t);
    local_f[1] = h * (f1 + 2 * f2) / 6.0 + elem.sigma * h * (q_prev[elem.number] + 2 * q_prev[elem.number + 1]) / (6 * delta_t);
}

// Сборка глоабальной СЛАУ (глобального вектора и матрицы)
void mfe::assembly_global_slae(mesh& mesh, const double t, const double delta_t)
{
    global_A->set_elem_to_null();
    global_f.clear();
    global_f.resize(q.size());

    if (method_to_solve == method::NEWTON)
    {
        A->set_elem_to_null();
        nl_f.clear();
        nl_f.resize(q.size());
    }

    for (uint32_t ielem = 0; ielem < mesh.get_fe_count(); ielem++)
    {
        // Собираем локальную матрицу
        build_local_matrix(mesh[ielem], delta_t);
        
        // Собираем локальынй вектор
        build_local_vector(mesh[ielem], t, delta_t);
        
        if (method_to_solve == method::NEWTON)
        {
            A->add_to_elem(0, ielem, local_A[0][0]);
            A->add_to_elem(0, ielem + 1, local_A[1][1]);
            A->add_to_elem(1, ielem, local_A[1][0]);
            A->add_to_elem(2, ielem, local_A[0][1]);

            nl_f[ielem] += local_f[0];
            nl_f[ielem + 1] += local_f[1];

            linearization_newton(mesh[ielem], delta_t);
        }

        // Добавляем в глобальную матрицу
        global_A->add_to_elem(0, ielem, local_A[0][0]);
        global_A->add_to_elem(0, ielem + 1, local_A[1][1]);
        global_A->add_to_elem(1, ielem, local_A[1][0]);
        global_A->add_to_elem(2, ielem, local_A[0][1]);

        // Добавляем в глобальный вектор
        global_f[ielem] += local_f[0];
        global_f[ielem + 1] += local_f[1];
    }
}

// Задание начального условия
void mfe::initial_condition(mesh& mesh)
{
    double t0 = mesh.get_time(0);
    q[0] = u(mesh[0].x, t0);
    for (uint32_t i = 0; i < mesh.get_fe_count(); i++)
        q[i + 1] = u(mesh[i].x_next, t0);
}

// Учет первого краевого условия
void mfe::first_boundary_condition(mesh& mesh, const double t, const double delta_t)
{
    // На диагонали ставим 1
    global_A->set_elem(0, 0, 1.0);
    global_A->set_elem(2, 0, 0.0);
    global_A->set_elem(0, mesh.get_fe_count(), 1.0);
    global_A->set_elem(1, mesh.get_fe_count() - 1, 0.0);

    // В вектор правой части ставим известные значения неизвестной функции
    global_f[0] = u(mesh[0].x, t);
    global_f[mesh.get_fe_count()] = u(mesh[mesh.get_fe_count() - 1].x_next, t);

    if (method_to_solve == method::NEWTON)
    {
        A->set_elem(0, 0, 1.0);
        A->set_elem(2, 0, 0.0);
        A->set_elem(0, mesh.get_fe_count(), 1.0);
        A->set_elem(1, mesh.get_fe_count() - 1, 0.0);

        nl_f[0] = u(mesh[0].x, t);
        nl_f[mesh.get_fe_count()] = u(mesh[mesh.get_fe_count() - 1].x_next, t);
    }
}

// Решить систему нелинейных уравнений
std::pair<uint32_t, double> mfe::solve(mesh& mesh, uint32_t max_iter, double eps, method method)
{
    method_to_solve = method;

    if (method_to_solve == method::NEWTON)
    {
        A = new matrix(mesh.get_fe_count() + 1, 3);
        nl_f.resize(mesh.get_fe_count() + 1);
    }

    initial_condition(mesh);

    solver slv;

    uint32_t iter = 0;

    // Решаем на каждом временном слое
    for (uint32_t i = 1; i < mesh.get_time_size(); i++)
    {
        double delta_t = mesh.get_time(i) - mesh.get_time(i - 1);

        for ( ; iter < max_iter; iter++)
        {
            assembly_global_slae(mesh, mesh.get_time(i), delta_t);
            first_boundary_condition(mesh, mesh.get_time(i), delta_t);

            if (iter != 0)
            {
                if (residual() < eps)
                    break;

                //double w = 0.0;
                //if (method_to_solve == method::SIMPLE_ITERATION)
                //    w = one_dimensional_search::minimize(*global_A, q, q_prev_min, global_f);
                //else
                //    w = one_dimensional_search::minimize(*A, q, q_prev_min, nl_f);

                //add_relax(w);

                //if (residual() < eps)
                //    break;
            }

            slv.solve_by_LU(*global_A, global_f, q);
            q_prev_min = q_prev;
            q_prev = q;
        }
    }
 
    auto res = std::make_pair(iter, error(mesh));

    //save(directory, res);

    return res;
}

// Добавка параметра релаксации
void mfe::add_relax(const double w)
{
    for (uint32_t i = 0; i < q.size(); i++)
        q[i] = w * q[i] + (1 - w) * q_prev_min[i];
}

// Расчет невязки
double mfe::residual()
{
    std::unique_ptr<std::vector<double>> Aq;
    if (method_to_solve == method::SIMPLE_ITERATION)
    {
        Aq = global_A->dot(q);
        for (uint32_t i = 0; i < q.size(); i++)
            (*Aq)[i] -= global_f[i];

        return norm(*Aq) / norm(global_f);
    }

    Aq = A->dot(q);
    for (uint32_t i = 0; i < q.size(); i++)
        (*Aq)[i] -= nl_f[i];

    return norm(*Aq) / norm(nl_f);
}

// Расчет относительной невязки
double mfe::error(mesh& mesh)
{
    double err = 0.0;
    double exact_norm = 0.0;

    exact.resize(q.size(), 0);

    exact[0] = u(mesh[0].x, mesh.get_time(mesh.get_time_size() - 1));
    for (uint32_t i = 0; i < mesh.get_fe_count(); i++)
        exact[i + 1] = u(mesh[i].x_next, mesh.get_time(mesh.get_time_size() - 1));

    for (uint32_t i = 0; i < q.size(); i++)
    {
        err += (exact[i] - q[i]) * (exact[i] - q[i]);
        exact_norm += exact[i] * exact[i];
    }

    return sqrt(err) / sqrt(exact_norm);
}

// Сохранить результат в файл
void mfe::save(std::string dir, std::pair<uint32_t, double>& result)
{
    std::ofstream res(dir + "result.txt");

    if (res.is_open())
    {
        res << std::left
            << std::setw(15) << "Iterations: " << result.first << std::endl
            << std::setw(15) << "Error: " << result.second << std::endl << std::endl;

        res << std::left << std::setw(15) << "q*"
            << std::setw(15) << "q"
            << std::setw(15) << "|q*-q|" << std::endl;
        res << "------------------------------------------" << std::endl;

        for (uint32_t i = 0; i < q.size(); i++)
        {
            res << std::left << std::setw(15) << exact[i]
                << std::setw(15) << q[i]
                << std::setw(15) << abs(exact[i] - q[i]) << std::endl;
        }

        exact.clear();
        res.close();
    }
}

//-------------------------------------------------------------------------------------------
//----------------------------- РЕАЛИЗАЦИЯ МЕТОДА НЬЮТОНА -----------------------------------
//
// Добавки к локальным матрицам и локальному вектору правой части
void mfe::linearization_newton(const finite_elem& elem, const double delta_t)
{
    double h = elem.x_next - elem.x;
    uint32_t ielem = elem.number;

    double dAii_dq1 =  1 / (2 * h) * (dlambda(q_prev[ielem + 1], q_prev[ielem], h, elem.x, 1)
                                    + dlambda(q_prev[ielem + 1], q_prev[ielem], h, elem.x_next, 1));
    double dAii_dq2 =  1 / (2 * h) * (dlambda(q_prev[ielem + 1], q_prev[ielem], h, elem.x, 2)
                                    + dlambda(q_prev[ielem + 1], q_prev[ielem], h, elem.x_next, 2));
    double dAij_dq1 = -1 / (2 * h) * (dlambda(q_prev[ielem + 1], q_prev[ielem], h, elem.x, 1)
                                    + dlambda(q_prev[ielem + 1], q_prev[ielem], h, elem.x_next, 1));
    double dAij_dq2 = -1 / (2 * h) * (dlambda(q_prev[ielem + 1], q_prev[ielem], h, elem.x, 2)
                                    + dlambda(q_prev[ielem + 1], q_prev[ielem], h, elem.x_next, 2));

    double dbi_dqi = elem.sigma * h / (3 * delta_t);
    double dbi_dqj = elem.sigma * h / (6 * delta_t);

    local_A[0][0] += dAii_dq1 * q_prev[ielem] + dAij_dq1 * q_prev[ielem + 1] - dbi_dqi;
    local_A[0][1] += dAii_dq2 * q_prev[ielem] + dAij_dq2 * q_prev[ielem + 1] - dbi_dqj;
    local_A[1][0] += dAij_dq1 * q_prev[ielem] + dAii_dq1 * q_prev[ielem + 1] - dbi_dqj;
    local_A[1][1] += dAij_dq2 * q_prev[ielem] + dAii_dq2 * q_prev[ielem + 1] - dbi_dqi;

    local_f[0] += q_prev[ielem] * (dAii_dq1 * q_prev[ielem] + dAii_dq2 * q_prev[ielem + 1])
                + q_prev[ielem + 1] * (dAij_dq1 * q_prev[ielem] + dAij_dq2 * q_prev[ielem + 1])
                - (dbi_dqi * q_prev[ielem] + dbi_dqj * q_prev[ielem + 1]);

    local_f[1] += q_prev[ielem] * (dAij_dq1 * q_prev[ielem] + dAij_dq2 * q_prev[ielem + 1])
                + q_prev[ielem + 1] * (dAii_dq1 * q_prev[ielem] + dAii_dq2 * q_prev[ielem + 1])
                - (dbi_dqj * q_prev[ielem] + dbi_dqi * q_prev[ielem + 1]);
}

// Производная от функции lambda по qi
double mfe::dlambda(double q2, double q1, double h, double x, uint32_t var)
{
    switch (var)
    {
    case 1:
        return (lambda(q2 / h, q1 / h + 1, x) - lambda(q2 / h, q1 / h, x)) / h;

    case 2:
        return (lambda(q2 / h + 1, q1 / h, x) - lambda(q2 / h, q1 / h, x)) / h;

    default:
        throw "Wrong variable";
        break;
    }
}