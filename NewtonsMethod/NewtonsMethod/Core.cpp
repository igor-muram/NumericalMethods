#include "Core.h"

const vector<function<double(double, double)>> Afuncs = {
	[](double x, double y) { return -3;		},
	[](double x, double y) { return 1;		},
	[](double x, double y) { return -6 * x; },
	[](double x, double y) { return 1;		}
};

const vector<function<double(double, double)>> Ffuncs = {
	[](double x, double y) { return y - 3 * x + 2; },
	[](double x, double y) { return y - 3 * x * x; }
};

void CalculateMatrix(double* A, double* x)
{
	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			A[i * N + j] = Afuncs[i * N + j](x[0], x[1]);
}

void CalculateF(double* F, double* x)
{
	for (int i = 0; i < N; i++)
		F[i] = Ffuncs[i](x[0], x[1]);
}

void Multiply(double* A, double* vec, double *res)
{
	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
		{
			int index = i * N + j;
			res[i] = A[index] * vec[j];
		}
}
