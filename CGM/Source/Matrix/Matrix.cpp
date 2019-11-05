#include "Matrix.h"					

void ReadMatrix(Matrix& A, int& Aaxiter, double& eps, int& choice)
{
	ifstream in("input/kuslau.txt");
	in >> A.N >> Aaxiter >> eps >> choice;
	in.close();

	A.IA = new int[A.N + 1];
	A.DI = new double[A.N];

	in.open("input/ig.txt");
	for (int i = 0; i < A.N + 1; i++)
		in >> A.IA[i];
	in.close();

	int size = A.IA[A.N];
	A.JA = new int[size];
	A.AL = new double[size];
	A.AU = new double[size];

	in.open("input/jg.txt");
	for (int i = 0; i < A.IA[A.N]; i++)
		in >> A.JA[i];
	in.close();

	in.open("input/ggl.txt");
	for (int i = 0; i < A.IA[A.N]; i++)
		in >> A.AL[i];
	in.close();

	in.open("input/ggu.txt");
	for (int i = 0; i < A.IA[A.N]; i++)
		in >> A.AU[i];
	in.close();

	in.open("input/di.txt");
	for (int i = 0; i < A.N; i++)
		in >> A.DI[i];
	in.close();
}

void ReadB(int N, double *b)
{
	ifstream in("input/pr.txt");
	for (int i = 0; i < N; i++)
		in >> b[i];
	in.close();
}

void ReadX0(int N, double *x0)
{
	ifstream in("x0.txt");

	int i = 0;
	for (; in >> x0[i] && i < N; i++);

	for (; i < N; i++)
		x0[i] = 0;
}

void Multiply(Matrix& A, double *vec, double *res)
{
	double *di = A.DI;
	double *al = A.AL;
	double *au = A.AU;
	int *ia = A.IA;
	int *ja = A.JA;
	int N = A.N;

	for (int i = 0; i < N; i++)
	{
		res[i] = vec[i] * di[i];
		for (int k = ia[i]; k < ia[i + 1]; k++)
		{
			int j = ja[k];
			res[i] += al[k] * vec[j];
			res[j] += au[k] * vec[i];
		}
	}
}

void MultiplyT(Matrix& A, double *vec, double *res)
{
	double *di = A.DI;
	double *al = A.AU;
	double *au = A.AL;
	int *ia = A.IA;
	int *ja = A.JA;
	int N = A.N;

	for (int i = 0; i < N; i++)
	{
		res[i] = vec[i] * di[i];
		for (int k = ia[i]; k < ia[i + 1]; k++)
		{
			int j = ja[k];
			res[i] += al[k] * vec[j];
			res[j] += au[k] * vec[i];
		}
	}
}

void LUFactorization(Matrix& A, Matrix& LU)
{
	LU.N = A.N;
	for (int i = 0; i < A.N + 1; i++)
		LU.IA[i] = A.IA[i];

	for (int i = 0; i < A.IA[A.N]; i++)
		LU.JA[i] = A.JA[i];

	for (int i = 0; i < A.N; i++)
	{
		double sumD = 0;
		int i0 = A.IA[i], i1 = A.IA[i + 1];
		int j = i - (i1 - i0);

		for (int k = i0; k < i1; k++, j++)
		{
			double sumL = 0, sumU = 0;

			// Calculate L[i][j], U[j][i]
			int j0 = A.IA[j], j1 = A.IA[j + 1];
			int size_i = k - i0, size_j = j1 - j0;
			int diff = size_i - size_j;
			int kl = i0, ku = j0;

			(diff < 0) ? ku -= diff : kl += diff;

			for (; kl < k; kl++, ku++)
			{
				sumL += A.AL[kl] * A.AU[ku];
				sumU += A.AU[kl] * A.AL[ku];
			}

			LU.AL[k] = A.AL[k] - sumL;
			LU.AU[k] = A.AU[k] - sumU;
			LU.AU[k] /= A.DI[j];

			// Calculate sum for DI[i]
			sumD += A.AL[k] * A.AU[k];
		}

		// Calculate DI[i]
		LU.DI[i] = A.DI[i] - sumD;
		sumD = 0;
	}
}

double DotProduct(int N, double *a, double *b)
{
	double sum = 0;
	for (int i = 0; i < N; i++)
		sum += a[i] * b[i];

	return sum;
}
