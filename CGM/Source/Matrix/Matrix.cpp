#include "Matrix.h"					

void ReadMatrix(Matrix& A, int& maxiter, double& eps, int& choice)
{
	std::ifstream in("input/kuslau.txt");
	in >> A.N >> maxiter >> eps >> choice;
	in.close();

	A.IA = new int[A.N + 1];
	A.DI = new double[A.N];

	in.open("input/ig.txt");
	for (int i = 0; i < A.N + 1; i++)
	{
		in >> A.IA[i];
		A.IA[i]--;
	}
	in.close();

	int size = A.IA[A.N];
	A.JA = new int[size];
	A.AL = new double[size];
	A.AU = new double[size];

	in.open("input/jg.txt");
	for (int i = 0; i < A.IA[A.N]; i++)
	{
		in >> A.JA[i];
		A.JA[i]--;
	}
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

void ReadB(int N, double* b)
{
	std::ifstream in("input/pr.txt");
	for (int i = 0; i < N; i++)
		in >> b[i];
	in.close();
}

void ReadX0(int N, double* x0)
{
	std::ifstream in("input/x0.txt");

	for (int i = 0; i < N; i++)
		in >> x0[i];

	in.close();
}

void Multiply(Matrix& A, double* vec, double* res)
{
	double* di = A.DI;
	double* al = A.AL;
	double* au = A.AU;
	int* ia = A.IA;
	int* ja = A.JA;
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

void MultiplyT(Matrix& A, double* vec, double* res)
{
	double* di = A.DI;
	double* al = A.AU;
	double* au = A.AL;
	int* ia = A.IA;
	int* ja = A.JA;
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

void MultiplyU(Matrix& A, double* vec, double* res)
{
	double* au = A.AU;
	int* ia = A.IA;
	int* ja = A.JA;
	int N = A.N;

	for (int i = 0; i < N; i++)
	{
		res[i] = vec[i];
		for (int k = ia[i]; k < ia[i + 1]; k++)
		{
			int j = ja[k];
			res[j] += au[k] * vec[i];
		}
	}
}

void LUFactorization(Matrix& A, Matrix& LU)
{
	LU.N = A.N;
	LU.IA = new int[LU.N + 1];

	for (int i = 0; i < A.N + 1; i++)
		LU.IA[i] = A.IA[i];

	LU.AL = new double[LU.IA[LU.N]];
	LU.AU = new double[LU.IA[LU.N]];
	LU.JA = new int[LU.IA[LU.N]];
	LU.DI = new double[LU.N];

	for (int i = 0; i < A.IA[A.N]; i++)
		LU.JA[i] = A.JA[i];

	for (int i = 0; i < A.N; i++)
	{
		double sumD = 0;
		int i0 = A.IA[i], i1 = A.IA[i + 1];

		for (int k = i0; k < i1; k++)
		{
			double sumL = 0, sumU = 0;
			int j = A.JA[k];

			// Calculate L[i][j], U[j][i]
			int j0 = A.IA[j], j1 = A.IA[j + 1];

			int kl = i0, ku = j0;

			for ( ; kl < i1 && ku < j1; )
			{
				int j_kl = A.JA[kl];
				int j_ku = A.JA[ku];

				if (j_kl == j_ku)
				{
					sumL += A.AL[kl] * A.AU[ku];
					sumU += A.AU[kl] * A.AL[ku];
					kl++;
					ku++;
				}
				if (j_kl > j_ku)
					ku++;
				if (j_kl < j_ku)
					kl++;
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

void Forward(Matrix& A, double* x, double* b, bool transposed)
{
	double* di = A.DI;
	double* al = A.AL;
	double* au = A.AU;
	int* ia = A.IA;
	int* ja = A.JA;
	int N = A.N;

	if (transposed)
	{
		for (int i = 0; i < N; i++)
			x[i] = b[i];

		for (int i = N - 1; i >= 0; i--)
		{
			int i0 = ia[i], i1 = ia[i + 1];
			x[i] /= di[i];
			for (int k = i0; k < i1; k++)
			{
				int j = ja[k];
				x[j] -= al[k] * x[i];

			}
		}
	}
	else
	{
		for (int i = 0; i < N; i++)
		{
			double sum = 0;
			int i0 = ia[i], i1 = ia[i + 1];
			for (int k = i0; k < i1; k++)
			{
				int j = ja[k];
				sum += al[k] * x[j];
			}
			x[i] = (b[i] - sum) / di[i];
		}
	}
}

void Backward(Matrix& A, double* x, double* b, bool transposed)
{
	double* di = A.DI;
	double* al = A.AL;
	double* au = A.AU;
	int* ia = A.IA;
	int* ja = A.JA;
	int N = A.N;

	if (transposed)
	{
		for (int i = 0; i < N; i++)
		{
			double sum = 0;
			int i0 = ia[i], i1 = ia[i + 1];
			for (int k = i0; k < i1; k++)
			{
				int j = ja[k];
				sum += au[k] * x[j];
			}
			x[i] = b[i] - sum;
		}
	}
	else
	{
		for (int i = 0; i < N; i++)
			x[i] = b[i];

		for (int i = N - 1; i >= 0; i--)
		{
			int i0 = ia[i], i1 = ia[i + 1];
			for (int k = i0; k < i1; k++)
			{
				int j = ja[k];
				x[j] -= au[k] * x[i];

			}
		}
	}
}

double DotProduct(int N, double* a, double* b)
{
	double sum = 0;
	for (int i = 0; i < N; i++)
		sum += a[i] * b[i];

	return sum;
}

Matrix HilbertMatrix(int size)
{
	Matrix A = { };
	int ALSize = size * (size - 1) / 2;

	A.N = size;
	A.DI = new double[size];
	A.AL = new double[ALSize];
	A.AU = new double [ALSize];
	A.IA = new int[size + 1];
	A.JA = new int[ALSize];

	for (int i = 0; i < size; i++)
		A.DI[i] = 1.0 / (2 * i + 1);

	A.IA[0] = 0;
	for (int i = 0; i < size; i++)
		A.IA[i + 1] = A.IA[i] + i;

	for (int i = 1; i < size; i++)
		for (int j = 0; j < i; j++)
		{
			A.JA[A.IA[i] + j] = j;
			A.AU[A.IA[i] + j] = A.AL[A.IA[i] + j] = 1.0 / (i + j + 1);
		}

	return A;
}
