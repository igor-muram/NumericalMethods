#include "Matrix.h"

void InitMatrix(Matrix& A, int N, int m)
{
	A.N = N;
	A.m = m;
	int mSize = A.N - A.m - 1;

	A.D.resize(A.N, 1);
	A.L1.resize(A.N - 1);
	A.U1.resize(A.N - 1);
	A.L2.resize(mSize);
	A.L2.resize(mSize);
}