#include "los.h"

int LOS(Matrix& A, double* x, double* f, double* r, double* z, double* p, double* Ax, int maxiter, double eps)
{
	// Ax = A * x
	Multiply(A, x, Ax);

	// r0 = f - Ax
	// z0 = r0
	for (int i = 0; i < A.N; i++)
	{
		r[i] = f[i] - Ax[i];
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

		Multiply(A, r, Ax);
		double bk = -DotProduct(A.N, p, Ax) / dotP;

		for (int i = 0; i < A.N; i++)
		{
			z[i] = r[i] + bk * z[i];
			p[i] = Ax[i] + bk * p[i];
		}

		diff = DotProduct(A.N, r, r);
	}

	return k;
}
