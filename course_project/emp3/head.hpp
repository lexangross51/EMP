#ifndef HEAD_HPP
#define HEAD_HPP

#include <iostream>
#include <vector>
#include <fstream>
#include <set>
#include "array"
#include <iomanip>
#include <functional>
#include <string>
#include <nlohmann/json.hpp>

const std::string directory = "C:\\Users\\lexan\\OneDrive\\������� ����\\����\\6 �������\\��������� �������������� ������\\course_project\\emp3\\";

typedef std::vector<double> dvector;

typedef std::function<double(double)> function1D;
typedef std::function<double(double, double)> function2D;
typedef std::function<double(double, double, double)> function2D_t;
typedef std::function<double(double, double, double, double, double, double)> function2D_t_r;

#endif