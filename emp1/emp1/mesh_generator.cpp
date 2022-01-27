#include "mesh_generator.h"

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
		file >> n_omega;
		areas.resize(n_omega);

		for (uint32_t i = 0; i < n_omega; i++)
		{
			uint32_t startx, endx, starty, endy;
			double lambda, gamma;
			file >> startx >> endx >> starty >> endy >> lambda >> gamma;

			point startp(X_lines[startx - 1], Y_lines[starty - 1]);
			point endp(X_lines[endx - 1], Y_lines[endy - 1]);

			areas[i].borders[0] = startp;
			areas[i].borders[1] = endp;
			areas[i].lambda = lambda;
			areas[i].gamma = gamma;
		}

		//==========================================
		file >> num_x >> num_y;
		file >> kx >> ky;

		file.close();
	}
	else throw "Can't open file";
}

// Генерируем массивы X и Y с учетом разбиения
void mesh_generator::generate_xy()
{
	double dx;
	// Если хотим неравномерную сетку
	if (kx != 1.0)
		dx = (X_lines[nx - 1] - X_lines[0]) * (1 - kx) / (1 - pow(kx, num_x));
	// Если хотим равномерную сетку
	else 
		dx = (X_lines[nx - 1] - X_lines[0]) / num_x;

	// Для оси y аналогично
	double dy;
	if (ky != 1.0)
		dy = (Y_lines[ny - 1] - Y_lines[0]) * (1 - ky) / (1 - pow(ky, num_y));
	else
		dy = (Y_lines[ny - 1] - Y_lines[0]) / num_y;

	// ================ Для оси X ================
	double coord = X_lines[0];

	X.push_back(coord);

	for (uint32_t i = 0, j = 1; i < num_x; i++)
	{
		coord += dx;
		dx *= kx;

		if (abs(coord - X_lines[j]) < 1e-14)
		{
			X.push_back(X_lines[j++]);
		}
		else if (coord > X_lines[j])
		{
			X.push_back(X_lines[j++]);
			X.push_back(coord);
		}
		else
			X.push_back(coord);
	}
	X_lines.clear();

	// ================ Для оси Y ================
	coord = Y_lines[0];

	Y.push_back(coord);

	for (uint32_t i = 0, j = 1; i < num_y; i++)
	{
		coord += dy;
		dy *= ky;

		if (abs(coord - Y_lines[j]) < 1e-14)
		{
			Y.push_back(Y_lines[j++]);
		}
		else if (coord > Y_lines[j])
		{
			Y.push_back(Y_lines[j++]);
			Y.push_back(coord);
		}
		else
			Y.push_back(coord);
	}
	Y_lines.clear();
}

// Строим сетку с учетом границ подобластей
void mesh_generator::build_mesh(std::vector<node>& mesh)
{
	input(directory);
	generate_xy();

	for (uint32_t i = 0; i < Y.size(); i++)
		for (uint32_t j = 0; j < X.size(); j++)
		{
			uint32_t area_num;

			point p(X[j], Y[i]);
			auto node_type = what_type(p);

			if (node_type == node::node_type::INTERNAL)
			{
				is_in_area(p, area_num);
				mesh.push_back(	node(p, node_type,
								areas[area_num].lambda,
								areas[area_num].gamma));
			}
			else 
				mesh.push_back(node(p, node_type, 0.0, 0.0));
		}

	X.clear();
	Y.clear();
	areas.clear();
	
	write_by_type(mesh, directory);
}

// Определяем тип узла: внутренний, граничный или фиктивный
node::node_type mesh_generator::what_type(point& p)
{
	uint32_t cnt = 0;

	uint32_t neighbor_cnt = 0;

	uint32_t ix = 0, iy = 0;

	// Определяем позицию точки в массивах X и Y
	for (uint32_t i = 0; i < X.size(); i++)
		if (p.x == X[i])
		{
			ix = i;
			break;
		}

	for (uint32_t i = 0; i < Y.size(); i++)
		if (p.y == Y[i])
		{
			iy = i;
			break;
		}

	// Проверим сначала, не является ли узел фиктивным
	// (т.е. не принадлежащим ни одной области)
	for (const auto& it : areas)
	{
		if (!(p.x >= it.borders[0].x && p.x <= it.borders[1].x &&
			p.y >= it.borders[0].y && p.y <= it.borders[1].y))
			cnt++;
	}
	if (cnt == areas.size())
		return node::node_type::FICTITIOUS;

	// Определим сколько соседей есть у точки
	// Возьмем узлы слева, справа, снизу и сверху от текущего
	uint32_t x_next = ix + 1;
	int x_prev = ix - 1;

	uint32_t y_next = iy + 1;
	int y_prev = iy - 1;

	// Если вышло так, что соседние узлы либо по X, либо по Y
	// вышли за пределы массивов X или Y, то узел был граничным
	if (x_prev < 0 || x_next > X.size() - 1 ||
		y_prev < 0 || y_next > Y.size() - 1)
		return node::node_type::BORDER;

	uint32_t c;

	if (is_in_area(point(X[x_prev], Y[iy]), c))
		neighbor_cnt++;
	if (is_in_area(point(X[x_next], Y[iy]), c))
		neighbor_cnt++;
	if (is_in_area(point(X[ix], Y[y_prev]), c))
		neighbor_cnt++;
	if (is_in_area(point(X[ix], Y[y_next]), c))
		neighbor_cnt++;

	// Если у узла 4 соседа, то он внутренний
	if (neighbor_cnt == 4)
		return node::node_type::INTERNAL;
	// Иначе граничный
	else
		return node::node_type::BORDER;
}

// Проверяет, находится ли точка в области
bool mesh_generator::is_in_area(const point& p, uint32_t &area_num)
{
	for (uint32_t i = 0; i < areas.size(); i++)
	{
		if (p.x >= areas[i].borders[0].x && p.x <= areas[i].borders[1].x &&
			p.y >= areas[i].borders[0].y && p.y <= areas[i].borders[1].y)
		{
			area_num = i;
			return true;
		}
	}
	return false;
}

// Записать узлы сетки в файл в зависимости от типа узла
void mesh_generator::write_by_type(const std::vector<node>& mesh, 
								   const std::string dir)
{
	std::ofstream i(dir + "internal.txt");
	std::ofstream b(dir + "border.txt");
	std::ofstream f(dir + "fictitious.txt");

	if (i.is_open() && b.is_open() && f.is_open())
	{
		for (const auto& it : mesh)
		{
			if (it.type == node::node_type::INTERNAL)
				i << it.p.x << " " << it.p.y << std::endl;
			else if (it.type == node::node_type::BORDER)
				b << it.p.x << " " << it.p.y << std::endl;
			else 
				f << it.p.x << " " << it.p.y << std::endl;
		}
		i.close();
		b.close();
		f.close();
	}
}