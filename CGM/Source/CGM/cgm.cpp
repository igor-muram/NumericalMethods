#include "cgm.h"

int CGM(Matrix& A, double *x, double *f, Memory& cache, int maxiter, double eps)
{
	double *Ax, *r, *z;
	cache.Allocate(&Ax);
	cache.Allocate(&r);
	cache.Allocate(&z);

	// Calculate r0, z0
	Multiply(A, x, Ax);
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
		// Calculate alpha
		Multiply(A, z, Ax);
		double a = dotR0 / DotProduct(A.N, Ax, z);

		// Calculate xk, rk
		for (int i = 0; i < A.N; i++)
		{
			x[i] += a * z[i];
			r[i] -= a * Ax[i];
		}

		// Calculate beta
		double dotR1 = DotProduct(A.N, r, r);
		double b = dotR1 / dotR0;

		// Calculate zk
		for (int i = 0; i < A.N; i++)
			z[i] = r[i] + b * z[i];

		// Calculate difference
		diff = dotR1 / dotF;
		dotR0 = dotR1;
	}

	cache.Free(&Ax);
	cache.Free(&r);
	cache.Free(&z);

	return k;
}

int CGM_diag(Matrix& A, double *x, double *f, Memory& cache, int maxiter, double eps)
{
	double *Ax, *r, *z, *L, *U, *temp;
	cache.Allocate(&Ax);
	cache.Allocate(&r);
	cache.Allocate(&z);
	cache.Allocate(&L);
	cache.Allocate(&U);
	cache.Allocate(&temp);

	for (int i = 0; i < A.N; i++)
		L[i] = U[i] = sqrt(fabs(A.DI[i]));

	// Calculate r0
	Multiply(A, x, Ax);
	for (int i = 0; i < A.N; i++)
		r[i] = (f[i] - Ax[i]) / (L[i] * L[i]);

	MultiplyT(A, r, temp);
	for (int i = 0; i < A.N; i++)
		r[i] = temp[i] / U[i];

	// Calculate z0
	for (int i = 0; i < A.N; i++)
		z[i] = r[i];

	//Calculate x0
	for (int i = 0; i < A.N; i++)
		x[i] *= U[i];

	double dotF = DotProduct(A.N, f, f);
	double dotR0 = DotProduct(A.N, r, r);
	double diff = dotR0 / dotF;

	int k = 0;
	for (; k < maxiter && diff >= eps * eps; k++)
	{
		// Calculate alpha
		for (int i = 0; i < A.N; i++)
			Ax[i] = z[i] / U[i];

		Multiply(A, Ax, temp);
		for (int i = 0; i < A.N; i++)
			Ax[i] = temp[i] / (L[i] * L[i]);

		MultiplyT(A, Ax, temp);
		for (int i = 0; i < A.N; i++)
			Ax[i] = temp[i] / U[i];

		double a = dotR0 / DotProduct(A.N, Ax, z);

		// Calculate xk, rk
		for (int i = 0; i < A.N; i++)
		{
			x[i] += a * z[i];
			r[i] -= a * Ax[i];
		}

		// Calculate beta
		double dotR1 = DotProduct(A.N, r, r);
		double b = dotR1 / dotR0;

		// Calculate zk
		for (int i = 0; i < A.N; i++)
			z[i] = r[i] + b * z[i];

		// Calculate difference
		diff = dotR1 / dotF;
		dotR0 = dotR1;
	}

	// Calculate final x
	for (int i = 0; i < A.N; i++)
		x[i] /= U[i];

	cache.Free(&Ax);
	cache.Free(&r);
	cache.Free(&z);
	cache.Free(&L);
	cache.Free(&U);
	cache.Free(&temp);

	return k;
}