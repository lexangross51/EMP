#include "mesh_generator.h"

// ������ ������ � ��������� �������
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

			areas[i].borders[0] = { {p1, p2},  border::bound_cond(bc1) };
			areas[i].borders[1] = { {p1, p3},  border::bound_cond(bc2) };
			areas[i].borders[2] = { {p3, p4},  border::bound_cond(bc3) };
			areas[i].borders[3] = { {p2, p4},  border::bound_cond(bc4) };
			areas[i].lambda = lambda;
			areas[i].gamma = gamma;
		}

		file.close();
	}
	else throw "Can't open file";
}

// ���������� ������� X � Y � ������ ���������
void mesh_generator::generate_xy(const grid_type type)
{
	double dx, dy;
	double coord;

	// ���� ����� ������ ���� �����������
	if (type == grid_type::UNIFORM)
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
			}
			else
				Y.push_back(coord);
		}

	}
	// ���� ����� ������ ���� �������������
	else
	{
		for (uint32_t i = 0; i < nx - 1; i++)
		{
			if (kr[0][i] != 1.0)
				dx = (X_lines[i + 1] - X_lines[i]) * (1 - kr[0][i]) / (1 - pow(kr[0][i], part[0][i]));
			else
				dx = (X_lines[nx - 1] - X_lines[0]) / part[0][0];
			
			coord = X_lines[i];
			
			for (; coord < X_lines[i + 1]; coord += dx)
			{
				X.push_back(coord);
				dx *= kr[0][i];
			}
		}
		X.push_back(X_lines[nx - 1]);

		for (uint32_t i = 0; i < part[1].size(); i++)
		{
			if (kr[1][i] != 1.0)
				dy = (Y_lines[i + 1] - Y_lines[i]) * (1 - kr[1][i]) / (1 - pow(kr[1][i], part[1][i]));
			else
				dy = (Y_lines[ny - 1] - Y_lines[0]) / part[1][0];

			coord = Y_lines[i];

			for (; coord < Y_lines[i + 1]; coord += dy)
			{
				Y.push_back(coord);
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

// ������ ����� � ������ ������ �����������
void mesh_generator::build_mesh(std::vector<node>& mesh, grid_type type)
{
	input(directory);
	generate_xy(type);

	for (uint32_t i = 0; i < Y.size(); i++)
		for (uint32_t j = 0; j < X.size(); j++)
		{
			point p(X[j], Y[i]);
			auto node_type = what_type(p, j, i);

			if (node_type != node::node_type::FICTITIOUS)
			{
				uint32_t area_num;
				is_in_area(p, area_num);

				if (node_type == node::node_type::INTERNAL)
					mesh.push_back(node(p, j, i, node_type,
						border::bound_cond::NONE,
						areas[area_num].lambda,
						areas[area_num].gamma));
				else
					mesh.push_back(node(p, j, i, node_type,
						areas[area_num].borders[what_border(p, area_num)].bc,
						areas[area_num].lambda,
						areas[area_num].gamma));
			}
			else
				mesh.push_back(node(p, 0, 0, node_type));
		}

	X.clear();
	Y.clear();
	areas.clear();
	
	write_by_type(mesh, directory);
}

// ���������� ��� ����: ����������, ��������� ��� ���������
node::node_type mesh_generator::what_type(const point& p, const uint32_t i, const uint32_t j)
{
	uint32_t cnt = 0;

	uint32_t neighbor_cnt = 0;

	// �������� �������, �� �������� �� ���� ���������
	for (const auto& it : areas)
	{
		if (!(p.x >= it.borders[0].limits[0].x && p.x <= it.borders[0].limits[1].x &&
			p.y >= it.borders[1].limits[0].y && p.y <= it.borders[1].limits[1].y))
			cnt++;
	}
	// ���� ����� ��������, ������� �� ���������� �����, ����� ������ ����� ��������
	// �� ���� - ���������
	if (cnt == areas.size())
		return node::node_type::FICTITIOUS;

	// ��������� ������� ������� ���� � �����
	// ������� ���� �����, ������, ����� � ������ �� ��������
	uint32_t x_next = i + 1;
	int x_prev = i - 1;

	uint32_t y_next = j + 1;
	int y_prev = j - 1;

	// ���� ����� ���, ��� �������� ���� ���� �� X, ���� �� Y
	// ����� �� ������� �������� X ��� Y, �� ���� ��� ���������
	if (x_prev < 0 || x_next > X.size() - 1 ||
		y_prev < 0 || y_next > Y.size() - 1)
		return node::node_type::BORDER;

	uint32_t c;
	if (is_in_area(point(X[x_prev], Y[j]), c))
		neighbor_cnt++;
	if (is_in_area(point(X[x_next], Y[j]), c))
		neighbor_cnt++;
	if (is_in_area(point(X[i], Y[y_prev]), c))
		neighbor_cnt++;
	if (is_in_area(point(X[i], Y[y_next]), c))
		neighbor_cnt++;

	// ���� � ���� 4 ������, �� �� ����������
	if (neighbor_cnt == 4)
		return node::node_type::INTERNAL;
	// ����� ���������
	else
		return node::node_type::BORDER;
}

// ���������� ����� �������, ������� ����������� �����
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
			return i;
	}
}

// ���������, ��������� �� ����� � ������� � ���������� ����� �������, 
// ���� ����� ����������� ��
bool mesh_generator::is_in_area(const point& p, uint32_t &area_num)
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

// �������� ���� ����� � ���� � ����������� �� ���� ����
void mesh_generator::write_by_type(const std::vector<node>& mesh, 
								   const std::string dir)
{
	std::ofstream i(dir + "internal.txt");
	std::ofstream b(dir + "border.txt");
	std::ofstream f(dir + "fictitious.txt");
	std::ofstream fi(dir + "first.txt");
	std::ofstream se(dir + "second.txt");
	std::ofstream th(dir + "third.txt");

	if (i.is_open() && b.is_open() && f.is_open() &&
		fi.is_open() && se.is_open() && th.is_open())
	{
		for (const auto& it : mesh)
		{
			if (it.type == node::node_type::INTERNAL)
				i << it.p.x << " " << it.p.y << std::endl;
			else if (it.type == node::node_type::BORDER)
			{
				b << it.p.x << " " << it.p.y << std::endl;

				if (it.bc == border::bound_cond::DIRICHLET)
					fi << it.p.x << " " << it.p.y << std::endl;
				else if (it.bc == border::bound_cond::NEUMANN)
					se << it.p.x << " " << it.p.y << std::endl;
				else if (it.bc == border::bound_cond::NEWTON)
					th << it.p.x << " " << it.p.y << std::endl;
			}
			else
				f << it.p.x << " " << it.p.y << std::endl;
		}
		i.close();
		b.close();
		f.close();
		fi.close();
		se.close();
		th.close();
	}
	else
		std::cerr << "Can't open files" << std::endl;
}