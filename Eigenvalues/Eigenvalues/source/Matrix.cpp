#include "Matrix.h"

void ReadMatrix(Matrix& A, int& maxiter, double& eps)
{
	std::ifstream in("input/size.txt");
	in >> A.N >> maxiter >> eps;
	in.close();

	A.IA = new int[A.N + 1];
	A.DI = new double[A.N];
	in.open("input/ia.txt");
	for (int i = 0; i < A.N + 1; i++)
		in >> A.IA[i];
	in.close();

	A.DI = new double[A.N];
	A.AL = new double[A.IA[A.N] + 1];
	A.AU = new double[A.IA[A.N] + 1];

	in.open("input/di.txt");
	for (int i = 0; i < A.N; i++)
		in >> A.DI[i];
	in.close();

	in.open("input/al.txt");
	for (int i = 0; i < A.IA[A.N] + 1; i++)
		in >> A.AL[i];
	in.close();

	in.open("input/au.txt");
	for (int i = 0; i < A.IA[A.N] + 1; i++)
		in >> A.AU[i];
	in.close();
}

void LUDecomposition(Matrix& A)
{
	for (int i = 0; i < A.N; i++)
	{
		double sumD = 0;
		int i0 = A.IA[i], i1 = A.IA[i + 1];
		int j = i - (i1 - i0);

		for (int k = i0; k < i1; k++, j++)
		{
			double sumL = 0, sumU = 0;

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

			A.AL[k] -= sumL;
			A.AU[k] -= sumU;
			A.AU[k] /= A.DI[j];

			sumD += A.AL[k] * A.AU[k];
		}

		A.DI[i] -= sumD;
		sumD = 0;
	}
}

void Solve(Matrix& A, double* b, double* x)
{
	double* y = new double[A.N];
	for (int i = 0; i < A.N; i++)
		y[i] = b[i];

	for (int i = 0; i < A.N; i++)
	{
		double sumL = 0;
		int j = i - 1;

		for (int k = A.IA[i + 1] - 1; k >= A.IA[i]; k--, j--)
			sumL += A.AL[k] * y[j];

		y[i] = (b[i] - sumL) / A.DI[i];
	}

	for (int i = A.N - 1; i >= 0; i--)
	{
		int j = i - 1;
		x[i] = y[i];
		for (int k = A.IA[i + 1] - 1; k >= A.IA[i]; k--, j--)
			y[j] -= A.AU[k] * x[i];
	}
}

void Multiply(Matrix& A, double* vector, double* res)
{
	for (int i = 0; i < A.N; i++)
		res[i] = vector[i] * A.DI[i];

	for (int i = 1; i < A.N; i++)
	{
		int i0 = A.IA[i];
		int i1 = A.IA[i + 1];
		int j = i - (i1 - i0);

		for (int k = i0; k < i1; k++, j++)
		{
			res[i] += vector[j] * A.AL[k];
			res[j] += vector[i] * A.AU[k];
		}
	}
}

double Norm(int N, double* vector)
{
	double sum = 0;
	for (int i = 0; i < N; i++)
		sum += vector[i] * vector[i];

	return sqrt(sum);
}

double MaxEigenValue(Matrix& A, double* x0, double* x1, int maxiter, double eps)
{
	int N = A.N;

	double lambda0 = 0, lambda1;
	double diff = 1;

	for (int k = 0; k < maxiter && diff > eps; k++)
	{
		Multiply(A, x0, x1);

		lambda1 = Norm(N, x1) / Norm(N, x0);
		diff = abs((lambda1 - lambda0) / lambda1);

		if (k % 10 == 0)
		{
			double norm = Norm(N, x1);
			for (int i = 0; i < N; i++)
				x1[i] /= norm;
		}

		lambda0 = lambda1;
		for (int i = 0; i < N; i++)
			x0[i] = x1[i];
	}

	return lambda0;
}

double MinEigenValue(Matrix& A, double* x0, double* x1, int maxiter, double eps)
{
	int N = A.N;

	double lambda0 = 0, lambda1;
	double diff = 1;

	LUDecomposition(A);
	for (int k = 0; k < maxiter && diff > eps; k++)
	{
		Solve(A, x0, x1);

		lambda1 = Norm(N, x1) / Norm(N, x0);
		diff = abs(lambda1 - lambda0);

		if (k % 10 == 0)
		{
			double norm = Norm(N, x1);
			for (int i = 0; i < N; i++)
				x1[i] /= norm;
		}

		lambda0 = lambda1;
		for (int i = 0; i < N; i++)
			x0[i] = x1[i];
	}

	return 1 / lambda0;
}