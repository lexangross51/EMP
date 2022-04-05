// ================ AREA.HPP ================
#pragma once 
#ifndef AREA_HPP
#define AREA_HPP

// ����� � ���������� ������������
struct point3D
{
    double x, y, z;

    point3D(double _x = 0, double _y = 0, double _z = 0) :
        x(_x), y(_y), z(_z) {};
};

// �������� ���� ��������� �������
struct area
{
    double x_start, x_end;
    double y_start, y_end;
    double z_start, z_end;

    double lambda, sigma, hi;
};

// ������� ����� (������������ ������� ��������� ��������)
struct omega
{
    point3D start_point;
    point3D end_point;
};

#endif