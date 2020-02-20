#pragma once

#include <vector>

using namespace std;

struct Matrix
{
	vector<double> D, L1, L2, U1, U2;
	int N, m;
};

void InitMatrix(Matrix& A, int N, int m);