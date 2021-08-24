#pragma once

#include "BCG.h"

struct AuxVectors
{
	vector<double> r, z, p, y, r0;
	vector<double> Ax, Ap, Az;
};

double DotProduct(vector<double>& a, vector<double>& b)
{
	int N = a.size();
	double sum = 0;
	for (int i = 0; i < N; i++)
		sum += a[i] * b[i];

	return sum;
}

int BCGSTABLU(Matrix& A, vector<double>& x, vector<double>& f, Matrix& LU, AuxVectors* aux, int maxiter, double eps)
{
	vector<double>& r0 = aux->r0;
	vector<double>& r = aux->r;
	vector<double>& z = aux->z;
	vector<double>& p = aux->p;
	vector<double>& y = aux->y;
	vector<double>& Ax = aux->Ax;
	vector<double>& Ap = aux->Ap;
	vector<double>& Az = aux->Az;
	int N = A.N;

	//  r0
	A.Multiply(x, Ax);

	for (int i = 0; i < N; i++)
		r0[i] = f[i] - Ax[i];

	LU.Forward(r0, r0);

	// rk = r0
	for (int i = 0; i < N; i++)
		r[i] = r0[i];

	// z0
	LU.Backward(z, r);

	// diff
	double diff = DotProduct(r, r);

	int k = 0;
	for (; k < maxiter && diff >= eps; k++)
	{
		// L(-1) * A * U(-1) * z(k - 1)
		LU.Backward(Ax, z);
		A.Multiply(Ax, Az);
		LU.Forward(Az, Az);
		
		// alpha(k)
		double dot = DotProduct(r, r0);
		double a = dot / DotProduct(Az, r0);

		// p(k)
		for (int i = 0; i < N; i++)
			p[i] = r[i] - a * Az[i];

		// L(-1) * A * U(-1) * p(k)
		LU.Backward(Ax, p);
		A.Multiply(Ax, Ap);
		LU.Forward(Ap, Ap);

		// gamma(k)
		double g = DotProduct(p, Ap) / DotProduct(Ap, Ap);

		// y(k)
		for (int i = 0; i < N; i++)
			y[i] += a * z[i] + g * p[i];

		// r(k)
		for (int i = 0; i < N; i++)
			r[i] = p[i] - g * Ap[i];

		// betta(k)
		double b = a * DotProduct(r, r0) / (g * dot);

		// z(k)
		for (int i = 0; i < N; i++)
			z[i] = r[i] + b * z[i] - b * g * Az[i];

		// diff
		diff = DotProduct(r, r);
	}

	LU.Backward(x, y);

	return k;
}

int BCG(Matrix& A, vector<double>& x, vector<double>& b)
{
	AuxVectors aux;

	aux.r.resize(A.N);
	aux.z.resize(A.N);
	aux.p.resize(A.N);
	aux.y.resize(A.N);
	aux.r0.resize(A.N);
	aux.Ax.resize(A.N);
	aux.Ap.resize(A.N);
	aux.Az.resize(A.N);

	Matrix LU;
	A.LUFactorization(LU);

	return BCGSTABLU(A, x, b, LU, &aux, 100000, 1.0e-8);
}
