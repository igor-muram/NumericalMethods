#pragma once
#include <functional>
#include <cmath>

#include "Classes.h"

double Gauss2(std::function<double(double, double)> f, Rect rect)
{
	double w = rect.x1 - rect.x0;
	double h = rect.z1 - rect.z0;

	double J = w * h / 4;

	double p1 = 1 / sqrt(3);
	double p2 = -1 / sqrt(3);

	double result =
		f(w * (p1 + 1) / 2 + rect.x0, h * (p1 + 1) / 2 + rect.z0) +
		f(w * (p1 + 1) / 2 + rect.x0, h * (p2 + 1) / 2 + rect.z0) +
		f(w * (p2 + 1) / 2 + rect.x0, h * (p1 + 1) / 2 + rect.z0) +
		f(w * (p2 + 1) / 2 + rect.x0, h * (p2 + 1) / 2 + rect.z0);

	return J * result;
}