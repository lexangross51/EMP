#pragma once
#include "head.h"

// �������: ������ ������� (2 �����) � ��� ��������
struct border
{
	enum class bound_cond
	{
		NONE,
		DIRICHLET,
		NEUMANN,
		NEWTON
	};

	point limits[2];
	bound_cond bc = bound_cond::NONE;
};

// ������������ �������: 4 ������� � ������������
struct area
{
	border borders[4];
	double lambda = 0.0;
	double gamma = 0.0;
};