// ================ HEAD.HPP ================
#ifndef HEAD_HPP
#define HEAD_HPP

#include <iostream>
#include <vector>
#include <array>
#include <set>
#include <fstream>
#include <iomanip>
#include <functional>
#include <string>
#include <chrono>
#include <nlohmann/json.hpp>

const std::string directory = "C:\\Users\\lexan\\OneDrive\\Рабочий стол\\НГТУ\\6 семестр\\Уравнения математической физики\\emp3\\emp3\\";

typedef std::vector<double> dvector;

typedef std::function<double(double, double, double)> function3D;
typedef std::function<double(double, double, double, double, double, double, double)> r_function3D;

#endif