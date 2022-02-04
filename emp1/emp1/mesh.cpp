#include "mesh.h" 

// Сохранить сетку
void mesh::save(const std::string dir)
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
		for (const auto& it : nodes)
		{
			if (it.type == 0)
				i << it.p.x << " " << it.p.y << std::endl;
			else if (it.type == 1 || it.type == 2 || it.type == 3 || it.type == 4)
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