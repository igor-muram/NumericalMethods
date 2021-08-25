#pragma once
#include "stdafx.h"

#pragma region Constants
const double w = 1.0e+00;
const double lambda = 1.0e+04;
const double sigma = 1.0e+04;
const double chi = 4.405e-11;
#pragma endregion

#pragma region Functions
inline double usin(double x, double y, double z)
{
	return x * x * x;
}

inline double ucos(double x, double y, double z)
{
	return x * x * x;
}

inline double divgrad_usin(double x, double y, double z)
{
	return 3 * 2 * x;
}

inline double divgrad_ucos(double x, double y, double z)
{
	return 3 * 2 * x;
}

inline double fsin(double x, double y, double z)
{
	return -divgrad_usin(x, y, z) * lambda - w * sigma * ucos(x, y, z) - w * w * chi * usin(x, y, z);
}

inline double fcos(double x, double y, double z)
{
	return -divgrad_ucos(x, y, z) * lambda + w * sigma * usin(x, y, z) - w * w * chi * ucos(x, y, z);
}
#pragma endregion

#pragma region Raw interval from file
struct RawInterval
{
	double begin, end;
	double q;
	int n;
};

inline std::istream& operator>>(std::istream& s, RawInterval& i)
{
	s >> i.begin >> i.end >> i.n >> i.q;
	return s;
}
#pragma endregion

#pragma region Interval prepeared for mesh building
struct Interval
{
	double begin, end;
	int n[2];
};

using IIterator = std::vector<Interval>::iterator;
#pragma endregion

#pragma region FiniteElement
struct FE
{
	double x[2];
	double y[2];
	double z[2];
	std::vector<int> verts;
};

using FEIterator = std::vector<FE>::iterator;
#pragma endregion

#pragma region Matrices
inline int a(int i)
{
	i++;
	return (i - 1) % 2;
}

inline int b(int i)
{
	i++;
	return ((i - 1) / 2) % 2;
}

inline int c(int i)
{
	i++;
	return (i - 1) / 4;
}

const std::vector<std::vector<double>> G1 = { { 1, -1 }, { -1, 1 } };
const std::vector<std::vector<double>> M1 = { { 1.0 / 3, 1.0 / 6 }, { 1.0 / 6, 1.0 / 3 } };

inline double M(int i, int j, double hx, double hy, double hz)
{
	return hx * hy * hz * M1[a(i)][a(j)] * M1[b(i)][b(j)] * M1[c(i)][c(j)];
}

inline double G(int i, int j, double hx, double hy, double hz)
{
	return
		hy * hz * G1[a(i)][a(j)] * M1[b(i)][b(j)] * M1[c(i)][c(j)] / hx +
		hx * hz * M1[a(i)][a(j)] * G1[b(i)][b(j)] * M1[c(i)][c(j)] / hy +
		hx * hy * M1[a(i)][a(j)] * M1[b(i)][b(j)] * G1[c(i)][c(j)] / hz;
}
#pragma endregion


struct Point
{
	double value;
	int i;
};
using PI = std::vector<Point>::iterator;