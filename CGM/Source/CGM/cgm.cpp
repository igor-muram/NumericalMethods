#include "cgm.h"

int CGM(Matrix& A, double* x0, double* f, double* r, double* z, int maxiter, double eps, double* x)
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

	double dotF = DotProduct(A.N, f, f);
	double diff = DotProduct(A.N, r, r) / dotF;

	int k = 0;
	for (; k < maxiter && diff >= eps * eps; k++)
	{
		// x0 = A * z(k - 1)
		Multiply(A, z, x0);
		double dotR0 = DotProduct(A.N, r, r);
		double ak = dotR0 / DotProduct(A.N, x0, z);

		for (int i = 0; i < A.N; i++)
		{
			x[i] += ak * z[i];
			r[i] -= ak * x0[i];
		}

		double dotR1 = DotProduct(A.N, r, r);
		double bk = dotR1 / dotR0;
		
		for (int i = 0; i < A.N; i++)
			z[i] = r[i] + bk * z[i];

		diff = dotR1 / dotF;
	}

	return k;
}
