#include "los.h"

int LOS(Matrix& A, double *x, double *f, Memory& cache, int maxiter, double eps)
{
	double *Ax, *r, *z, *p;
	cache.Allocate(&Ax);
	cache.Allocate(&r);
	cache.Allocate(&z);
	cache.Allocate(&p);

	// Calculate r0, z0
	Multiply(A, x, Ax);
	for (int i = 0; i < A.N; i++)
	{
		r[i] = f[i] - Ax[i];
		z[i] = r[i];
	}

	// Calculate p0
	Multiply(A, z, p);

	double diff = DotProduct(A.N, r, r);

	int k = 0;
	for (; k < maxiter && diff >= eps; k++)
	{
		// Calculate alpha
		double dotP = DotProduct(A.N, p, p);
		double a = DotProduct(A.N, p, r) / dotP;

		// Calculate xk, rk
		for (int i = 0; i < A.N; i++)
		{
			x[i] += a * z[i];
			r[i] -= a * p[i];
		}

		// Calculate beta
		Multiply(A, r, Ax);
		double b = -DotProduct(A.N, p, Ax) / dotP;

		// Calculate zk, pk
		for (int i = 0; i < A.N; i++)
		{
			z[i] = r[i] + b * z[i];
			p[i] = Ax[i] + b * p[i];
		}

		// Calculate difference
		diff = DotProduct(A.N, r, r);
	}

	return k;
}

int LOS_diag(Matrix& A, double *x, double *f, Memory& cache, int maxiter, double eps)
{
	double *Ax, *r, *z, *p, *L, *U, *temp;
	cache.Allocate(&Ax);
	cache.Allocate(&r);
	cache.Allocate(&z);
	cache.Allocate(&p);
	cache.Allocate(&L);
	cache.Allocate(&U);
	cache.Allocate(&temp);

	for (int i = 0; i < A.N; i++)
		L[i] = U[i] = sqrt(fabs(A.DI[i]));

	// Calculate r0, z0
	Multiply(A, x, Ax);
	for (int i = 0; i < A.N; i++)
	{
		r[i] = (f[i] - Ax[i]) / L[i];
		z[i] = r[i] / U[i];
	}

	// Calculate p0
	Multiply(A, z, p);
	for (int i = 0; i < A.N; i++)
		p[i] /= L[i];

	double diff = DotProduct(A.N, r, r);

	int k = 0;
	for (; k < maxiter && diff >= eps; k++)
	{
		// Calculate alpha
		double dotP = DotProduct(A.N, p, p);
		double a = DotProduct(A.N, p, r) / dotP;

		// Calculate xk, rk
		for (int i = 0; i < A.N; i++)
		{
			x[i] += a * z[i];
			r[i] -= a * p[i];
		}

		// Calculate beta
		for (int i = 0; i < A.N; i++)
			temp[i] = r[i] / U[i];

		Multiply(A, temp, Ax);
		for (int i = 0; i < A.N; i++)
			Ax[i] /= L[i];

		double b = -DotProduct(A.N, p, Ax) / dotP;

		// Calculate zk, pk
		for (int i = 0; i < A.N; i++)
		{
			z[i] = r[i] / U[i] + b * z[i];
			p[i] = Ax[i] + b * p[i];
		}

		// Calculate difference
		diff = DotProduct(A.N, r, r);
	}

	cache.Free(&Ax);
	cache.Free(&r);
	cache.Free(&z);
	cache.Free(&p);
	cache.Free(&L);
	cache.Free(&U);
	cache.Free(&temp);

	return k;
}