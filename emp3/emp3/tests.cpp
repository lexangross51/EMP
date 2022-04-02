#include "tests.hpp"

function3D us = [](double x, double y, double z) { return x * x * x + y * y * y + z * z * z; };
function3D uc = [](double x, double y, double z) { return 2 * x * x * x - y * y * y + 3 * z * z * z; };

r_function3D fs = [](double x, double y, double z, double lambda, double sigma, double hi, double w)
{
	return -lambda * 6 * (x + y + z) - w * sigma * uc(x, y, z) - w * w * hi * us(x, y, z);
};

r_function3D fc = [](double x, double y, double z, double lambda, double sigma, double hi, double w)
{
	return -lambda * 6 * (2 * x - y + 3 * z) + w * sigma * us(x, y, z) - w * w * hi * uc(x, y, z);
};


void tests::calc_exact(space_grid& grid)
{
	exact.resize(2 * grid.get_nodes_count());

	for (uint32_t i = 0; i < grid.get_nodes_count(); i++)
	{
		point3D point = grid.get_point(i);

		exact[2 * i] = us(point.x, point.y, point.z);
		exact[2 * i + 1] = uc(point.x, point.y, point.z);
	}
}


void tests::omega_tests()
{
	space_grid_generator sgg;
	space_grid* sg = nullptr;

	sgg.build_mesh(sg, us, uc);

	calc_exact(*sg);

	double omegas[] = { 1e-4, 1e-3, 1e-2, 1e-1, 1, 10, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9 };

	result res1, res2, res3;

	std::ofstream wres(directory + "tests\\" + "omega_tests.txt");
	std::ofstream wresdif(directory + "tests\\" + "omega_tests_diffs.txt");

	wres << std::left << std::setw(10) << "omega"
		<< std::setw(10) << "|LU_time" << std::setw(10) << "|LOS_time" << std::setw(14) << "|BCGSTAB_time"
		<< std::setw(14) << "|LOS_iters" << std::setw(16) << "|BCGSTAB_iters" <<
		std::setw(16) << "|LOS_residual" << std::setw(16) << "|BSG_residual" << std::endl;
	wres << "-------------------------------------------------------------------------------------------------------" << std::endl;

	for (const auto& w : omegas)
	{
		mfe mfe(*sg);
		mfe.set_w(w);
		mfe.assembly_global_matrix_and_vector(fs, fc);

		res1 = mfe.solve(method::LU);
		print_diff(mfe.q, wresdif, "LU");

		res2 = mfe.solve(method::LOS_LU);
		print_diff(mfe.q, wresdif, "LOS");

		res3 = mfe.solve(method::BCGSTAB_LU);
		print_diff(mfe.q, wresdif, "BCGSTAB");

		wres << std::left << std::setw(10) << w
			<< "|" << std::setw(9) << res1.time.count() << "|" << std::setw(9) << res2.time.count() << "|" << std::setw(13) << res3.time.count()
			<< "|" << std::setw(13) << res2.iters << "|" << std::setw(15) << res3.iters <<
			"|" << std::setw(15) << res2.residual << "|" << std::setw(15) << res3.residual << std::endl;

		mfe.~mfe();
	}
}

void tests::lambda_tests()
{
	space_grid_generator sgg;
	space_grid* sg = nullptr;

	sgg.build_mesh(sg, us, uc);

	calc_exact(*sg);

	double lambdas[] = { 1e2, 1e3, 1e4, 1e5, 8e5 };

	result res1, res2, res3;

	std::ofstream wres(directory + "tests\\" + "lambda_tests.txt");
	std::ofstream wresdif(directory + "tests\\" + "lambda_tests_diffs.txt");

	wres << std::left << std::setw(10) << "lambda"
		<< std::setw(10) << "|LU_time" << std::setw(10) << "|LOS_time" << std::setw(14) << "|BCGSTAB_time"
		<< std::setw(14) << "|LOS_iters" << std::setw(16) << "|BCGSTAB_iters" <<
		std::setw(16) << "|LOS_residual" << std::setw(16) << "|BSG_residual" << std::endl;
	wres << "-------------------------------------------------------------------------------------------------------" << std::endl;

	for (const auto& lambda : lambdas)
	{
		sg->set_lambda(lambda);

		mfe mfe(*sg);
		mfe.set_w(10);
		mfe.assembly_global_matrix_and_vector(fs, fc);

		res1 = mfe.solve(method::LU);
		print_diff(mfe.q, wresdif, "LU");

		res2 = mfe.solve(method::LOS_LU);
		print_diff(mfe.q, wresdif, "LOS");

		res3 = mfe.solve(method::BCGSTAB_LU);
		print_diff(mfe.q, wresdif, "BCGSTAB");

		wres << std::left << std::setw(10) << lambda
			<< "|" << std::setw(9) << res1.time.count() << "|" << std::setw(9) << res2.time.count() << "|" << std::setw(13) << res3.time.count()
			<< "|" << std::setw(13) << res2.iters << "|" << std::setw(15) << res3.iters <<
			"|" << std::setw(15) << res2.residual << "|" << std::setw(15) << res3.residual << std::endl;

		mfe.~mfe();
	}
}

void tests::sigma_tests()
{
	space_grid_generator sgg;
	space_grid* sg = nullptr;

	sgg.build_mesh(sg, us, uc);

	calc_exact(*sg);

	double sigmas[] = { 0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8 };

	result res1, res2, res3;

	std::ofstream wres(directory + "tests\\" + "sigma_tests.txt");
	std::ofstream wresdif(directory + "tests\\" + "sigma_tests_diffs.txt");

	wres << std::left << std::setw(10) << "sigma"
		<< std::setw(10) << "|LU_time" << std::setw(10) << "|LOS_time" << std::setw(14) << "|BCGSTAB_time"
		<< std::setw(14) << "|LOS_iters" << std::setw(16) << "|BCGSTAB_iters" <<
		std::setw(16) << "|LOS_residual" << std::setw(16) << "|BSG_residual" << std::endl;
	wres << "-------------------------------------------------------------------------------------------------------" << std::endl;

	for (const auto& sigma : sigmas)
	{
		sg->set_sigma(sigma);

		mfe mfe(*sg);
		mfe.set_w(10);
		mfe.assembly_global_matrix_and_vector(fs, fc);

		res1 = mfe.solve(method::LU);
		print_diff(mfe.q, wresdif, "LU");

		res2 = mfe.solve(method::LOS_LU);
		print_diff(mfe.q, wresdif, "LOS");

		res3 = mfe.solve(method::BCGSTAB_LU);
		print_diff(mfe.q, wresdif, "BCGSTAB");

		wres << std::left << std::setw(10) << sigma
			<< "|" << std::setw(9) << res1.time.count() << "|" << std::setw(9) << res2.time.count() << "|" << std::setw(13) << res3.time.count()
			<< "|" << std::setw(13) << res2.iters << "|" << std::setw(15) << res3.iters <<
			"|" << std::setw(15) << res2.residual << "|" << std::setw(15) << res3.residual << std::endl;


		mfe.~mfe();
	}
}

void tests::hi_tests()
{
	space_grid_generator sgg;
	space_grid* sg = nullptr;

	sgg.build_mesh(sg, us, uc);

	calc_exact(*sg);

	double his[] = { 8.81e-12, 1e-12, 1e-11, 1e-10 };

	result res1, res2, res3;

	std::ofstream wres(directory + "tests\\" + "hi_tests.txt");
	std::ofstream wresdif(directory + "tests\\" + "hi_tests_diffs.txt");

	wres << std::left << std::setw(10) << "hi"
		<< std::setw(10) << "|LU_time" << std::setw(10) << "|LOS_time" << std::setw(14) << "|BCGSTAB_time"
		<< std::setw(14) << "|LOS_iters" << std::setw(16) << "|BCGSTAB_iters" <<
		std::setw(16) << "|LOS_residual" << std::setw(16) << "|BSG_residual" << std::endl;
	wres << "-------------------------------------------------------------------------------------------------------" << std::endl;

	for (const auto& hi : his)
	{
		sg->set_sigma(hi);

		mfe mfe(*sg);
		mfe.set_w(10);
		mfe.assembly_global_matrix_and_vector(fs, fc);

		res1 = mfe.solve(method::LU);
		print_diff(mfe.q, wresdif, "LU");

		res2 = mfe.solve(method::LOS_LU);
		print_diff(mfe.q, wresdif, "LOS");

		res3 = mfe.solve(method::BCGSTAB_LU);
		print_diff(mfe.q, wresdif, "BCGSTAB");

		wres << std::left << std::setw(10) << hi
			<< "|" << std::setw(9) << res1.time.count() << "|" << std::setw(9) << res2.time.count() << "|" << std::setw(13) << res3.time.count()
			<< "|" << std::setw(13) << res2.iters << "|" << std::setw(15) << res3.iters <<
			"|" << std::setw(15) << res2.residual << "|" << std::setw(15) << res3.residual << std::endl;

		mfe.~mfe();
	}
}


void tests::print_diff(dvector& solution, std::ofstream& file, std::string method)
{
	file << "-----------------------" << method << "-----------------------" << std::endl;
	file << std::left << std::setw(14) << "u*" << std::setw(14) << "u" << std::setw(14) << "|u*-u|" << std::endl;

	for (uint32_t i = 0; i < exact.size(); i++)
		file << std::left << std::setw(14) << exact[i] 
		<< std::setw(14) << solution[i] 
		<< std::setw(14) << abs(exact[i] - solution[i]) << std::endl;

	auto dif = solution - exact;

	file << norm(dif) / norm(exact) << std::endl;

	file << std::endl << std::endl;
}


void tests::run()
{
	omega_tests();

	//lambda_tests();

	//sigma_tests();

	//hi_tests();
}
