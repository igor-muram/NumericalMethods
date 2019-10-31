#include "los.h"

int LOS(Matrix& A, double* x0, double* f, double* r, double* z, double* p, int maxiter, double eps, double* x)
{
	// Copy x0 to x
	for (int i = 0; i < A.N; i++)
		x[i] = x0[i];

	// x0 = A * x0
	Multiply(A, x0, x0);

	// r0 = f - x0
	// z0 = r0
	for (int i = 0; i < A.N; i++)
	{
		r[i] = f[i] - x0[i];
		z[i] = r[i];
	}

	// p0 = A * z0
	Multiply(A, z, p);

	double diff = DotProduct(A.N, r, r);
	int k = 0;
	for (; k < maxiter && diff >= eps; k++)
	{
		double dotP = DotProduct(A.N, p, p);
		double ak = DotProduct(A.N, p, r) / dotP;

		for (int i = 0; i < A.N; i++)
		{
			x[i] += ak * z[i];
			r[i] -= ak * p[i];
		}

		Multiply(A, r, x0);
		double bk = -DotProduct(A.N, p, x0) / dotP;

		for (int i = 0; i < A.N; i++)
		{
			z[i] = r[i] + bk * z[i];
			p[i] = x0[i] + bk * p[i];
		}

		diff = DotProduct(A.N, r, r);
	}

	return k;
}
