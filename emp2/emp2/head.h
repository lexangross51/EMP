#pragma once
#include <iostream>
#include <vector>
#include <functional>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>

typedef std::function<double(double, double)> function;
typedef std::function<double(double, double, double)> function_f;
typedef std::function<double(double, double, double)> function_lambda;

const std::string directory = "C:\\Users\\lexan\\OneDrive\\������� ����\\����\\6 �������\\��������� �������������� ������\\emp2\\emp2\\";
const std::string tests_directory = "C:\\Users\\lexan\\OneDrive\\������� ����\\����\\6 �������\\��������� �������������� ������\\emp2\\emp2\\tests\\";

// ��������� ����� �������
inline double norm(const std::vector<double>& v)
{
	double scalar = 0.0;

	for (const auto& it : v)
		scalar += it * it;

	return sqrt(scalar);
}