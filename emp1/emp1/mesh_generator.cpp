#include "mesh_generator.h"

// Читаем данные о расчетной области
void mesh_generator::input(const std::string dir)
{
    n_omega = nx = ny = 0;

    std::ifstream file(dir + "areas.txt");

    if (file.is_open())
    {
        //==========================================
        file >> nx >> ny;
        X_lines.resize(nx);
        Y_lines.resize(ny);

        for (uint32_t i = 0; i < nx; i++)
            file >> X_lines[i];

        for (uint32_t i = 0; i < ny; i++)
            file >> Y_lines[i];

        //==========================================
        part.resize(2);
        part[0].resize(nx - 1);
        part[1].resize(ny - 1);

        kr.resize(2);
        kr[0].resize(nx - 1);
        kr[1].resize(ny - 1);

        for (uint32_t i = 0; i < nx - 1; i++)
            file >> part[0][i];
        for (uint32_t i = 0; i < ny - 1; i++)
            file >> part[1][i];

        for (uint32_t i = 0; i < nx - 1; i++)
            file >> kr[0][i];
        for (uint32_t i = 0; i < ny - 1; i++)
            file >> kr[1][i];

        //==========================================
        file >> n_omega;
        areas.resize(n_omega);

        for (uint32_t i = 0; i < n_omega; i++)
        {
            uint32_t startx, endx, starty, endy;
            double lambda, gamma;
            uint16_t bc1, bc2, bc3, bc4;

            file >> startx >> endx >> starty >> endy
                 >> lambda >> gamma
                 >> bc1 >> bc2 >> bc3 >> bc4;

            point p1(X_lines[startx - 1], Y_lines[starty - 1]);
            point p2(X_lines[endx - 1], Y_lines[starty - 1]);
            point p3(X_lines[startx - 1], Y_lines[endy - 1]);
            point p4(X_lines[endx - 1], Y_lines[endy - 1]);

            areas[i].borders[0] = { {p1, p2}, border::bound_cond(bc1) };
            areas[i].borders[1] = { {p1, p3}, border::bound_cond(bc2) };
            areas[i].borders[2] = { {p3, p4}, border::bound_cond(bc3) };
            areas[i].borders[3] = { {p2, p4}, border::bound_cond(bc4) };
            areas[i].lambda = lambda;
            areas[i].gamma = gamma;
        }

        file.close();
    }
    else throw "Can't open file";
}

// Генерируем массивы X и Y с учетом разбиения
void mesh_generator::generate_xy(mesh::mesh_type& type)
{
    double dx, dy;
    double coord;

    // Если сетка должна быть равномерной
    if (type == mesh::mesh_type::UNIFORM)
    {
        dx = (X_lines[nx - 1] - X_lines[0]) / part[0][0];
        dy = (Y_lines[ny - 1] - Y_lines[0]) / part[1][0];

        //==============================================
        coord = X_lines[0];
        X.push_back(coord);

        for (uint32_t i = 0, j = 1; i < part[0][0]; i++)
        {
            coord += dx;

            if (abs(coord - X_lines[j]) < 1e-14)
            {
                X.push_back(X_lines[j++]);
            }
            else if (coord > X_lines[j])
            {
                X.push_back(X_lines[j++]);
                X.push_back(coord);
                type = mesh::mesh_type::NONUNIFORM;
            }
            else
            {
                X.push_back(coord);
            }
        }

        //==============================================
        coord = Y_lines[0];
        Y.push_back(coord);

        for (uint32_t i = 0, j = 1; i < part[1][0]; i++)
        {
            coord += dy;

            if (abs(coord - Y_lines[j]) < 1e-14)
            {
                Y.push_back(Y_lines[j++]);
            }
            else if (coord > Y_lines[j])
            {
                Y.push_back(Y_lines[j++]);
                Y.push_back(coord);
                type = mesh::mesh_type::NONUNIFORM;
            }
            else
                Y.push_back(coord);
        }

    }
    // Если сетка должна быть неравномерной
    else
    {
        for (uint32_t i = 0; i < nx - 1; i++)
        {
            if (kr[0][i] != 1.0)
                dx = (X_lines[i + 1] - X_lines[i]) * (1 - kr[0][i]) / (1 - pow(kr[0][i], part[0][i]));
            else
                dx = (X_lines[i + 1] - X_lines[i]) / part[0][i];
            
            coord = X_lines[i];
            for (uint32_t j = 0; j < part[0][i]; j++)
            {
                X.push_back(coord);
                coord += dx;
                dx *= kr[0][i];
            }
        }
        X.push_back(X_lines[nx - 1]);

        for (uint32_t i = 0; i < ny - 1; i++)
        {
            if (kr[1][i] != 1.0)
                dy = (Y_lines[i + 1] - Y_lines[i]) * (1 - kr[1][i]) / (1 - pow(kr[1][i], part[1][i]));
            else
                dy = (Y_lines[i + 1] - Y_lines[i]) / part[1][i];

            coord = Y_lines[i];
            for (uint32_t j = 0; j < part[1][i]; j++)
            {
                Y.push_back(coord);
                coord += dy;
                dy *= kr[1][i];
            }
        }
        Y.push_back(Y_lines[ny - 1]);
    }

    X_lines.clear();
    Y_lines.clear();
    kr.clear();
    part.clear();
}

// Строим сетку с учетом границ подобластей
void mesh_generator::build_mesh(mesh& mesh, mesh::mesh_type type)
{
    input(directory);
    generate_xy(type);

    mesh.set_type(type);
    mesh.set_width(X.size() - 1);
    mesh.set_height(Y.size() - 1);

    for (uint32_t i = 0; i < Y.size(); i++)
        for (uint32_t j = 0; j < X.size(); j++)
        {
            point p(X[j], Y[i]);
            int node_type = what_type(p, j, i);

            if (node_type != -1)
            {
                uint32_t area_num;
                is_in_area(p, area_num);

                if (node_type == 0)
                    mesh.add_node(node(p, j, i, node_type,
                        border::bound_cond::NONE,
                        areas[area_num].lambda,
                        areas[area_num].gamma));
                else
                    mesh.add_node(node(p, j, i, node_type,
                        areas[area_num].borders[what_border(p, area_num) - 1].bc,
                        areas[area_num].lambda,
                        areas[area_num].gamma));
            }
            else
                mesh.add_node(node(p, j, i, node_type));
        }

    X.clear();
    Y.clear();
    areas.clear();

    mesh.save();
}

// Определяем тип узла: внутренний, граничный или фиктивный
int mesh_generator::what_type(const point& p, const uint32_t i, const uint32_t j)
{
    uint32_t cnt = 0;

    uint32_t neighbor_cnt = 0;

    uint32_t c;

    // Проверим сначала, не является ли узел фиктивным
    for (const auto& it : areas)
        if (!is_in_area(p, c))
            cnt++;

    // Если число областей, которым не принадлжит точка, равно общему числу областей
    // то узле - фиктивный
    if (cnt == areas.size())
        return -1;

    // Определим сколько соседей есть у точки
    // Возьмем узлы слева, справа, снизу и сверху от текущего
    uint32_t x_next = i + 1;
    int x_prev = i - 1;

    uint32_t y_next = j + 1;
    int y_prev = j - 1;

    // Если вышло так, что соседние узлы либо по X, либо по Y
    // вышли за пределы массивов X или Y, то узел был граничным
    if (x_prev < 0 || x_next > X.size() - 1 ||
        y_prev < 0 || y_next > Y.size() - 1)
    {
        is_in_area(p, c);
        return what_border(p, c);
    }

    if (is_in_area(point(X[x_prev], Y[j]), c))
        neighbor_cnt++;
    if (is_in_area(point(X[x_next], Y[j]), c))
        neighbor_cnt++;
    if (is_in_area(point(X[i], Y[y_prev]), c))
        neighbor_cnt++;
    if (is_in_area(point(X[i], Y[y_next]), c))
        neighbor_cnt++;

    // Если у узла 4 соседа, то он внутренний
    if (neighbor_cnt == 4)
        return 0;
    // Иначе граничный
    else
    {
        is_in_area(p, c);
        return what_border(p, c);
    }
}

// Возвращает номер границы, которой принадлежит точка
uint32_t mesh_generator::what_border(const point& p, const uint32_t area_num)
{
    for (uint32_t i = 0; i < 4; i++)
    {
        if (p.x == areas[area_num].borders[i].limits[0].x &&
            p.y >= areas[area_num].borders[i].limits[0].y &&
            p.y <= areas[area_num].borders[i].limits[1].y ||
            p.y == areas[area_num].borders[i].limits[0].y &&
            p.x >= areas[area_num].borders[i].limits[0].x &&
            p.x <= areas[area_num].borders[i].limits[1].x)
            return (i + 1);
    }
}

// Проверяет, находится ли точка в области и возвращает номер области, 
// если точка принадлежит ей
bool mesh_generator::is_in_area(const point& p, uint32_t& area_num)
{
    for (uint32_t i = 0; i < areas.size(); i++)
    {
        if (p.x >= areas[i].borders[0].limits[0].x &&
            p.x <= areas[i].borders[0].limits[1].x &&
            p.y >= areas[i].borders[1].limits[0].y &&
            p.y <= areas[i].borders[1].limits[1].y)
        {
            area_num = i;
            return true;
        }
    }
    return false;
}