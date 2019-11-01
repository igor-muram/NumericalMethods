#include "cgm.h"

int CGM(Matrix& A, double* x, double* f, double* r, double* z, double* Ax, int maxiter, double eps)
{
	// Ax = A * x0
	Multiply(A, x, Ax);

	// r0 = f - Ax
	// z0 = r0
	for (int i = 0; i < A.N; i++)
	{
		r[i] = f[i] - Ax[i];
		z[i] = r[i];
	}

	double dotF = DotProduct(A.N, f, f);
	double dotR0 = DotProduct(A.N, r, r);
	double diff = dotR0 / dotF;

	int k = 0;
	for (; k < maxiter && diff >= eps * eps; k++)
	{
		Multiply(A, z, Ax);
		double ak = dotR0 / DotProduct(A.N, Ax, z);

		for (int i = 0; i < A.N; i++)
		{
			x[i] += ak * z[i];
			r[i] -= ak * Ax[i];
		}

		double dotR1 = DotProduct(A.N, r, r);
		double bk = dotR1 / dotR0;
		
		for (int i = 0; i < A.N; i++)
			z[i] = r[i] + bk * z[i];

		diff = dotR1 / dotF;

		dotR0 = dotR1;
	}

	return k;
}
