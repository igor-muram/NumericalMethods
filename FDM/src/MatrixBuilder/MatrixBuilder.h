#pragma once

#include <vector>
#include <functional>

#include <SLAE/SLAE.h>
#include "../MeshBuilder/MeshBuilder.h"

using namespace std;

const double gamma = 2.0;

const function<double(double, double)> f = [](double x, double y) { return -2 + gamma * x * x; };

const vector<function<double(double, double)>> borderFuncs =
{
	[](double x, double y) { return x * x; },
	[](double x, double y) { return 1.0; },
	[](double x, double y) { return 4.0; },
	[](double x, double y) { return x * x; },
	[](double x, double y) { return x * x; },
	[](double x, double y) { return 4.0; },
	[](double x, double y) { return 9.0; },
	[](double x, double y) { return x * x; }
};

void BuildMatrix(
	SLAE::Matrix& A,
	vector<double>& b,
	vector<vector<int>>& areas,
	vector<Interval>& intervalsX,
	vector<Interval>& intervalsY,
	vector<BoundaryCondition>& conds,
	vector<double>& x,
	vector<double>& y,
	int kx, int ky);

void BoundaryConditions(
	SLAE::Matrix& A,
	vector<double>& b,
	vector<double>& x,
	vector<double>& y,
	vector<BoundaryCondition>& conds);

int IntervalNo(vector<Interval>& intervals, int index);
bool IsOnBorder(vector<BoundaryCondition>& conds, int ix, int iy);
