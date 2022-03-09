#pragma once
#ifndef MFE_H
#define MFE_H

#include "mesh.h"
#include "matrix.h"
#include "solver.h"

class mfe
{
	friend class solver;

public:
	enum class method
	{
		SIMPLE_ITERATION,
		NEWTON
	};

	mfe(mesh& mesh);

	void set_functions(const function& _u, const function_f& _f, const function_lambda& _lambda);

	void assembly_global_slae(mesh& mesh, const double t, const double delta_t);

	void initial_condition(mesh& mesh);
	void first_boundary_condition(mesh& mesh, const double t, const double delta_t);

	std::pair<uint32_t, double> solve(mesh& mesh, uint32_t max_iter, double eps, method method);

private:
	std::vector<std::vector<double>> local_A;	// ��������� �������
	std::vector<double> local_f;				// ��������� ������
	
	matrix* global_A;				// ���������� �������
	std::vector<double> global_f;	// ���������� ������

	matrix* A;						// ���������� ������� ��� ������� �������, ���� 
									// ���������� ����� �������
	std::vector<double> nl_f;

	std::vector<double> q;			// ������ ����� �� ������� ��������� ����
	std::vector<double> q_prev;		// ������ ����� �� ���������� ��������� ����

	std::vector<double> exact;		// ������ ��������

	function u;				// ����������� �������
	function_f f;			// ������� ������ �����
	function_lambda lambda;	// ������� lambda(du/dx)

	method method_to_solve;

	double dlambda(double q2, double q1, double h, double x, uint32_t var);

	void build_local_matrix(const finite_elem& elem, const double delta_t);
	void build_local_vector(const finite_elem& elem, const double t, const double delta_t);

	void add_relax(const double w);

	double residual();
	double error(mesh& mesh);

	void save(std::string dir, std::pair<uint32_t, double>& res);

	void linearization_newton(const finite_elem& elem, const double delta_t);
};

#endif 