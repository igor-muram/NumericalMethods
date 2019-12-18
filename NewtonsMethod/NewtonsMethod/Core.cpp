#include "Core.h"

const vector<vector<function<double(double, double)>>> Afuncs = {
	{
		[](double x, double y) { return 2 * (x + 1); },
		[](double x, double y) { return 2 * (y - 1); },
	},
	{
		[](double x, double y) { return 2 * (x - 3); },
		[](double x, double y) { return 2 * (y - 1); }
	},
	{
		[](double x, double y) { return 1; },
		[](double x, double y) { return -1; }
	}
};

const vector<function<double(double, double)>> Ffuncs = {
	[](double x, double y) { return  (x + 1) * (x + 1) + (y - 1) * (y - 1) - 4; },
	[](double x, double y) { return (x - 3) * (x - 3) + (y - 1) * (y - 1) - 4; },
	[](double x, double y) { return (x - y); }
};

void CalculateMatrix(double** A, double* x)
{
	for (int i = 0; i < M; i++)
		for (int j = 0; j < N; j++)
			A[i][j] = Afuncs[i][j](x[0], x[1]);
}

void CalculateF(double* F, double* x)
{
	for (int i = 0; i < M; i++)
		F[i] = -Ffuncs[i](x[0], x[1]);
}

void CalculateF(double* F, double* x, double* dx, double beta)
{
	for (int i = 0; i < M; i++)
		F[i] = Ffuncs[i](x[0] + beta * dx[0], x[1] + beta * dx[1]);
}

void ProcessSLAE(double** A, double* F)
{
	for (int i = 0; i < M; i++)
		for (int j = i + 1; j < M; j++)
			if (abs(F[i]) < abs(F[j]))
			{
				swap(F[i], F[j]);
				for (int k = 0; k < N; k++)
					swap(A[i][k], A[j][k]);
			}
}

void Multiply(double** A, double* vec, double* res)
{
	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			res[i] = A[i][j] * vec[j];
}

void Solve(double** A, double* x, double* F)
{
	for (int k = 0; k < N - 1; k++)
	{
		for (int i = k + 1; i < N; i++)
		{
			double t = A[i][k] / A[k][k];
			F[i] -= t * F[k];
			for (int j = k + 1; j < N; j++)
				A[i][j] -= t * A[k][j];
		}
	}

	x[N - 1] = F[N - 1] / A[N - 1][N - 1];
	for (int k = N - 2; k >= 0; k--)
	{
		double sum = 0;
		for (int j = k + 1; j < N; j++)
			sum += A[k][j] * x[j];

		x[k] = (F[k] - sum) / A[k][k];
	}
}

double Norm(double* vec)
{
	double sum = 0;
	for (int i = 0; i < M; i++)
		sum += vec[i] * vec[i];

	return sqrt(sum);
}