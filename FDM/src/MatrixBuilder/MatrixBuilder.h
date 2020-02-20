#pragma once

#include <vector>

#include "../MeshBuilder/MeshBuilder.h"
#include "Matrix.h"

using namespace std;

const double gamma = 2.0;

void BuildMatrix(
	Matrix& A,
	vector<vector<int>>& areas,
	vector<Interval>& intervalsX,
	vector<Interval>& intervalsY,
	vector<double>& x,
	vector<double>& y,
	vector<double>& hx,
	vector<double>& hy,
	int kx, int ky);

int IntervalNo(vector<Interval>& intervals, int index);