#pragma once

#include <vector>
#include <functional>

struct FiniteElement
{
	double begin, end;
	std::vector<int> v;
};

using FEIterator = std::vector<FiniteElement>::iterator;

std::vector<std::vector<double>> M = {
	{ 2.0 / 15, 1.0 / 15, -1.0 / 30 },
	{ 1.0 / 15, 8.0 / 15, 1.0 / 15 },
	{ -1.0 / 30, 1.0 / 15, 2.0 / 15 }
};

std::vector<std::vector<double>> G1 = {
	{ 37.0 / 30, -22.0 / 15, 7.0 / 30 },
	{ -22.0 / 15, 8.0 / 5, -2.0 / 15 },
	{ 7.0 / 30, -2.0 / 15, -1.0 / 10 }
};

std::vector<std::vector<double>> G2 = {
	{ 6.0 / 5, -16.0 / 15, -4.0 / 30 },
	{ -16.0 / 15, 32.0 / 15, -16.0 / 15 },
	{ -4.0 / 30, -16.0 / 15, 6.0 / 5 }
};

std::vector<std::vector<double>> G3 = {
	{ -1.0 / 10, -2.0 / 15, 7.0 / 30 },
	{ -2.0 / 15, 8.0 / 5, -22.0 / 15 },
	{ 7.0 / 30, -22.0 / 15, 37.0 / 30 }
};

std::vector<std::function<double(double, double)>> basis =
{
	[](double x, double h) { return 2 * (x - h) * (x - h / 2) / (h * h); },
	[](double x, double h) { return 4 * x * (h - x) / (h * h); },
	[](double x, double h) { return 2 * x * (x - h / 2) / (h * h); }
};

std::vector<std::function<double(double, double)>> basisDirs =
{
	[](double x, double h) { return 4 * x / (h * h) - 3 / h; },
	[](double x, double h) { return 4 / h - 8 * x / (h * h); },
	[](double x, double h) { return 4 * x / (h * h) - 1 / h; }
};

const double sigma = 0.0;

double Lamdba(double x, double du)
{
	return 1.0;
}

double F(double x, double t)
{
	return 0.0;
}

const double PI = 3.14159265359;