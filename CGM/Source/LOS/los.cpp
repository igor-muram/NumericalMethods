#include "los.h"

int LOS(Matrix& A, double* x, double* f, AuxVectors& aux, int maxiter, double eps, double& lastdiff)
{
	double* Ax = aux.Ax;
	double* r = aux.r;
	double* z = aux.z;
	double* p = aux.p;

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

	lastdiff = diff;
	return k;
}

int LOS_diag1(Matrix& A, double* x, double* f, AuxVectors& aux, int maxiter, double eps, double& lastdiff)
{
	double* Ax = aux.Ax;
	double* r = aux.r;
	double* z = aux.z;
	double* p = aux.p;
	double* LU = aux.LU;
	double* temp = aux.temp;
	int N = A.N;

	for (int i = 0; i < N; i++)
		LU[i] = sqrt(A.DI[i]);

	// Calculate r0, z0
	Multiply(A, x, Ax);
	for (int i = 0; i < N; i++)
	{
		r[i] = (f[i] - Ax[i]) / LU[i];
		z[i] = r[i] / LU[i];
	}

	// Calculate p0
	Multiply(A, z, p);
	for (int i = 0; i < N; i++)
		p[i] /= LU[i];

	double diff = DotProduct(N, r, r);

	int k = 0;
	for (; k < maxiter && diff >= eps; k++)
	{
		// Calculate alpha
		double dotP = DotProduct(N, p, p);
		double a = DotProduct(N, p, r) / dotP;

		// Calculate xk, rk
		for (int i = 0; i < N; i++)
		{
			x[i] += a * z[i];
			r[i] -= a * p[i];
		}

		// Calculate beta
		for (int i = 0; i < N; i++)
			temp[i] = r[i] / LU[i];

		Multiply(A, temp, Ax);
		for (int i = 0; i < N; i++)
			Ax[i] /= LU[i];

		double b = -DotProduct(N, p, Ax) / dotP;

		// Calculate zk, pk
		for (int i = 0; i < N; i++)
		{
			z[i] = r[i] / LU[i] + b * z[i];
			p[i] = Ax[i] + b * p[i];
		}

		// Calculate difference
		diff = DotProduct(A.N, r, r);
	}

	lastdiff = diff;
	return k;
}

int LOS_diag2(Matrix& A, double* x, double* f, AuxVectors& aux, int maxiter, double eps, double& lastdiff)
{
	double* Ax = aux.Ax;
	double* r = aux.r;
	double* z = aux.z;
	double* p = aux.p;
	double* LU = aux.LU;
	int N = A.N;

	for (int i = 0; i < N; i++)
		LU[i] = A.DI[i];

	// Calculate r0, z0
	Multiply(A, x, Ax);
	for (int i = 0; i < N; i++)
		z[i] = r[i] = (f[i] - Ax[i]) / LU[i];

	// Calculate p0
	Multiply(A, z, p);
	for (int i = 0; i < N; i++)
		p[i] /= LU[i];

	double diff = DotProduct(N, r, r);

	int k = 0;
	for (; k < maxiter && diff >= eps; k++)
	{
		// Calculate alpha
		double dotP = DotProduct(N, p, p);
		double a = DotProduct(N, p, r) / dotP;

		// Calculate xk, rk
		for (int i = 0; i < N; i++)
		{
			x[i] += a * z[i];
			r[i] -= a * p[i];
		}

		// Calculate beta
		Multiply(A, r, Ax);
		for (int i = 0; i < N; i++)
			Ax[i] /= LU[i];

		double b = -DotProduct(N, p, Ax) / dotP;

		// Calculate zk, pk
		for (int i = 0; i < N; i++)
		{
			z[i] = r[i] + b * z[i];
			p[i] = Ax[i] + b * p[i];
		}

		// Calculate difference
		diff = DotProduct(N, r, r);
	}

	lastdiff = diff;
	return k;
}

int LOS_diag3(Matrix& A, double* x, double* f, AuxVectors& aux, int maxiter, double eps, double& lastdiff)
{
	double* Ax = aux.Ax;
	double* r = aux.r;
	double* z = aux.z;
	double* p = aux.p;
	double* LU = aux.LU;
	double* temp = aux.temp;
	int N = A.N;

	for (int i = 0; i < N; i++)
		LU[i] = A.DI[i];

	// Calculate r0, z0
	Multiply(A, x, Ax);
	for (int i = 0; i < N; i++)
	{
		r[i] = f[i] - Ax[i];
		z[i] = r[i] / LU[i];
	}

	// Calculate p0
	Multiply(A, z, p);

	double diff = DotProduct(N, r, r);

	int k = 0;
	for (; k < maxiter && diff >= eps; k++)
	{
		// Calculate alpha
		double dotP = DotProduct(N, p, p);
		double a = DotProduct(N, p, r) / dotP;

		// Calculate xk, rk
		for (int i = 0; i < N; i++)
		{
			x[i] += a * z[i];
			r[i] -= a * p[i];
		}

		// Calculate beta
		for (int i = 0; i < N; i++)
			temp[i] = r[i] / LU[i];

		Multiply(A, temp, Ax);
		double b = -DotProduct(N, p, Ax) / dotP;

		// Calculate zk, pk
		for (int i = 0; i < N; i++)
		{
			z[i] = r[i] / LU[i] + b * z[i];
			p[i] = Ax[i] + b * p[i];
		}

		// Calculate difference
		diff = DotProduct(N, r, r);
	}

	lastdiff = diff;
	return k;
}

int LOS_LU(Matrix& A, double* x, double* f, Matrix& LU, AuxVectors& aux, int maxiter, double eps, double& lastdiff)
{
	double* Ax = aux.Ax;
	double* r = aux.r;
	double* z = aux.z;
	double* p = aux.p;
	double* temp = aux.temp;
	int N = A.N;

	// Calculate r0
	Multiply(A, x, Ax);
	for (int i = 0; i < N; i++)
		r[i] = f[i] - Ax[i];
	Forward(LU, r, r, false);

	//Calculate z0
	Backward(LU, z, r, false);

	// Calculate p0
	Multiply(A, z, p);
	Forward(LU, p, p, false);

	double diff = DotProduct(N, r, r);

	int k = 0;
	for (; k < maxiter && diff >= eps; k++)
	{
		// Calculate alpha
		double dotP = DotProduct(N, p, p);
		double a = DotProduct(N, p, r) / dotP;

		// Calculate xk, rk
		for (int i = 0; i < N; i++)
		{
			x[i] += a * z[i];
			r[i] -= a * p[i];
		}

		// Calculate beta
		Backward(LU, Ax, r, false);
		Multiply(A, Ax, temp);
		Forward(LU, Ax, temp, false);
		double b = -DotProduct(N, p, Ax) / dotP;

		// Calculate zk, pk
		Backward(LU, temp, r, false);
		for (int i = 0; i < N; i++)
		{
			z[i] = temp[i] + b * z[i];
			p[i] = Ax[i] + b * p[i];
		}

		// Calculate difference
		diff = DotProduct(N, r, r);
	}

	lastdiff = diff;
	return k;
}