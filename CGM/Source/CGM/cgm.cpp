#include "cgm.h"

int CGM(Matrix& A, double* x, double* f, AuxVectors& aux, int maxiter, double eps)
{
	double* Ax = aux.Ax;
	double* r = aux.r;
	double* z = aux.z;
	int N = A.N;

	// Calculate r0, z0
	Multiply(A, x, Ax);
	for (int i = 0; i < N; i++)
	{
		r[i] = f[i] - Ax[i];
		z[i] = r[i];
	}

	double dotF = DotProduct(N, f, f);
	double dotR0 = DotProduct(N, r, r);
	double diff = dotR0 / dotF;

	int k = 0;
	for (; k < maxiter && diff >= eps * eps; k++)
	{
		// Calculate alpha
		Multiply(A, z, Ax);
		double a = dotR0 / DotProduct(N, Ax, z);

		// Calculate xk, rk
		for (int i = 0; i < N; i++)
		{
			x[i] += a * z[i];
			r[i] -= a * Ax[i];
		}

		// Calculate beta
		double dotR1 = DotProduct(N, r, r);
		double b = dotR1 / dotR0;

		// Calculate zk
		for (int i = 0; i < N; i++)
			z[i] = r[i] + b * z[i];

		// Calculate difference
		diff = dotR1 / dotF;
		dotR0 = dotR1;
	}

	return k;
}

int CGM_diag1(Matrix& A, double* x, double* f, AuxVectors& aux, int maxiter, double eps)
{
	double* Ax = aux.Ax;
	double* r = aux.r;
	double* z = aux.z;
	double* LU = aux.LU;
	double* temp = aux.temp;
	int N = A.N;

	for (int i = 0; i < N; i++)
		LU[i] = sqrt(A.DI[i]);

	// Calculate r0
	Multiply(A, x, Ax);
	for (int i = 0; i < N; i++)
		r[i] = (f[i] - Ax[i]) / (LU[i] * LU[i]);

	MultiplyT(A, r, temp);
	for (int i = 0; i < N; i++)
		r[i] = temp[i] / LU[i];

	// Calculate z0
	for (int i = 0; i < N; i++)
		z[i] = r[i];

	//Calculate x0
	for (int i = 0; i < N; i++)
		x[i] *= LU[i];

	double dotF = DotProduct(N, f, f);
	double dotR0 = DotProduct(N, r, r);
	double diff = dotR0 / dotF;

	int k = 0;
	for (; k < maxiter && diff >= eps * eps; k++)
	{
		// Calculate alpha
		for (int i = 0; i < N; i++)
			Ax[i] = z[i] / LU[i];

		Multiply(A, Ax, temp);
		for (int i = 0; i < N; i++)
			Ax[i] = temp[i] / (LU[i] * LU[i]);

		MultiplyT(A, Ax, temp);
		for (int i = 0; i < N; i++)
			Ax[i] = temp[i] / LU[i];

		double a = dotR0 / DotProduct(N, Ax, z);

		// Calculate xk, rk
		for (int i = 0; i < N; i++)
		{
			x[i] += a * z[i];
			r[i] -= a * Ax[i];
		}

		// Calculate beta
		double dotR1 = DotProduct(N, r, r);
		double b = dotR1 / dotR0;

		// Calculate zk
		for (int i = 0; i < N; i++)
			z[i] = r[i] + b * z[i];

		// Calculate difference
		diff = dotR1 / dotF;
		dotR0 = dotR1;
	}

	// Calculate final x
	for (int i = 0; i < N; i++)
		x[i] /= LU[i];

	return k;
}

int CGM_diag2(Matrix& A, double* x, double* f, AuxVectors& aux, int maxiter, double eps)
{
	double* Ax = aux.Ax;
	double* r = aux.r;
	double* z = aux.z;
	double* LU = aux.LU;
	double* temp = aux.temp;
	int N = A.N;

	for (int i = 0; i < N; i++)
		LU[i] = A.DI[i];

	// Calculate r0
	Multiply(A, x, Ax);
	for (int i = 0; i < N; i++)
		temp[i] = (f[i] - Ax[i]) / (LU[i] * LU[i]);

	MultiplyT(A, temp, r);

	// Calculate z0
	for (int i = 0; i < N; i++)
		z[i] = r[i];

	double dotF = DotProduct(N, f, f);
	double dotR0 = DotProduct(N, r, r);
	double diff = dotR0 / dotF;

	int k = 0;
	for (; k < maxiter && diff >= eps * eps; k++)
	{
		// Calculate alpha
		Multiply(A, z, temp);
		for (int i = 0; i < N; i++)
			temp[i] = temp[i] / (LU[i] * LU[i]);

		MultiplyT(A, temp, Ax);
		double a = dotR0 / DotProduct(N, Ax, z);

		// Calculate xk, rk
		for (int i = 0; i < N; i++)
		{
			x[i] += a * z[i];
			r[i] -= a * Ax[i];
		}

		// Calculate beta
		double dotR1 = DotProduct(N, r, r);
		double b = dotR1 / dotR0;

		// Calculate zk
		for (int i = 0; i < N; i++)
			z[i] = r[i] + b * z[i];

		// Calculate difference
		diff = dotR1 / dotF;
		dotR0 = dotR1;
	}

	return k;
}

int CGM_diag3(Matrix& A, double* x, double* f, AuxVectors& aux, int maxiter, double eps)
{
	double* Ax = aux.Ax;
	double* r = aux.r;
	double* z = aux.z;
	double* LU = aux.LU;
	double* temp = aux.temp;
	int N = A.N;

	for (int i = 0; i < N; i++)
		LU[i] = A.DI[i];

	// Calculate r0
	Multiply(A, x, Ax);
	for (int i = 0; i < N; i++)
		r[i] = f[i] - Ax[i];

	MultiplyT(A, r, temp);
	for (int i = 0; i < N; i++)
		r[i] = temp[i] / LU[i];

	// Calculate z0
	for (int i = 0; i < N; i++)
		z[i] = r[i];

	//Calculate x0
	for (int i = 0; i < N; i++)
		x[i] *= LU[i];

	double dotF = DotProduct(N, f, f);
	double dotR0 = DotProduct(N, r, r);
	double diff = dotR0 / dotF;

	int k = 0;
	for (; k < maxiter && diff >= eps * eps; k++)
	{
		// Calculate alpha
		for (int i = 0; i < N; i++)
			Ax[i] = z[i] / LU[i];

		Multiply(A, Ax, temp);
		MultiplyT(A, temp, Ax);
		for (int i = 0; i < N; i++)
			Ax[i] /= LU[i];

		double a = dotR0 / DotProduct(N, Ax, z);

		// Calculate xk, rk
		for (int i = 0; i < N; i++)
		{
			x[i] += a * z[i];
			r[i] -= a * Ax[i];
		}

		// Calculate beta
		double dotR1 = DotProduct(N, r, r);
		double b = dotR1 / dotR0;

		// Calculate zk
		for (int i = 0; i < N; i++)
			z[i] = r[i] + b * z[i];

		// Calculate difference
		diff = dotR1 / dotF;
		dotR0 = dotR1;
	}

	// Calculate final x
	for (int i = 0; i < N; i++)
		x[i] /= LU[i];

	return k;
}