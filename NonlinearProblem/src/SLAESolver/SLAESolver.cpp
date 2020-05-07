#include "SLAESolver.h"

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

			// Вычисление элементов L[i][j] и U[j][i]
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

			// Накопление суммы для DI[i]
			sumD += A.AL[k] * A.AU[k];
		}

		// Вычисление DI[i]
		A.DI[i] -= sumD;
		sumD = 0;
	}

}

void Solve(Matrix& A, std::vector<double>& x, std::vector<double>& b)
{

	std::vector<double>& y = b;

	// Прямой обход
	for (int i = 0; i < A.N; i++)
	{
		double sumL = 0;
		int j = i - 1;

		for (int k = A.IA[i + 1] - 1; k >= A.IA[i]; k--, j--)
			sumL += A.AL[k] * y[j];

		y[i] = (b[i] - sumL) / A.DI[i];
	}

	// Обратный обход
	for (int i = A.N - 1; i >= 0; i--)
	{
		int j = i - 1;
		x[i] = y[i];
		for (int k = A.IA[i + 1] - 1; k >= A.IA[i]; k--, j--)
			y[j] -= A.AU[k] * x[i];
	}
}

void Multiply(Matrix& A, std::vector<double>& vector, std::vector<double>& res)
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