#include "Core.h"

const vector<function<double(double, double)>> F1 = {
	[](double x, double y) { return  y - x * x; },
	[](double x, double y) { return y - x - 2; },
};

const vector<vector<function<double(double, double)>>> A1 = {
	{
		[](double x, double y) { return -2 * x; },
		[](double x, double y) { return 1; }
	},
	{
		[](double x, double y) { return -1; },
		[](double x, double y) { return 1; }
	}
};

const vector<function<double(double, double)>> F2 = {
	[](double x, double y) { return  (x - 2) * (x - 2) + (y - 1) * (y - 1) - 4; },
	[](double x, double y) { return  (x + 2) * (x + 2) + (y - 1) * (y - 1) - 4; }
};

const vector<vector<function<double(double, double)>>> A2 = {
	{
		[](double x, double y) { return 2 * (x - 2); },
		[](double x, double y) { return 2 * (y - 1); },
	},
	{
		[](double x, double y) { return 2 * (x + 2); },
		[](double x, double y) { return 2 * (y - 1); }
	}
};

const vector<function<double(double, double)>> F3 = {
	[](double x, double y) { return  (x - 2) * (x - 2) + (y - 1) * (y - 1) - 4; },
	[](double x, double y) { return  (x + 2) * (x + 2) + (y - 1) * (y - 1) - 4; },
	[](double x, double y) { return   y - x - 1; }
};

const vector<vector<function<double(double, double)>>> A3 = {
	{
		[](double x, double y) { return 2 * (x - 2); },
		[](double x, double y) { return 2 * (y - 1); },
	},
	{
		[](double x, double y) { return 2 * (x + 2); },
		[](double x, double y) { return 2 * (y - 1); }
	},
	{
		[](double x, double y) { return -1; },
		[](double x, double y) { return 1; }
	}
};

const vector<vector<function<double(double, double, double)>>> Afuncs3 = {
	{
		[](double x, double y, double z) { return 4; },
		[](double x, double y, double z) { return 0; },
		[](double x, double y, double z) { return 4 * z; },
	},
	{
		[](double x, double y, double z) { return 0; },
		[](double x, double y, double z) { return -2; },
		[](double x, double y, double z) { return 3; }
	}
};

const vector<function<double(double, double, double)>> Ffuncs3 = {
	[](double x, double y, double z) { return  4 * x + 2 * z * z; },
	[](double x, double y, double z) { return  3 * z - 2 * y; }
};

void CalculateMatrix(double** A, double* x)
{
	for (int i = 0; i < M; i++)
		for (int j = 0; j < N; j++)
			A[i][j] = A1[i][j](x[0], x[1]);
}

void CalculateF(double* F, double* x)
{
	for (int i = 0; i < M; i++)
		F[i] = -F1[i](x[0], x[1]);
}

void CalculateF(double* F, double* x, double* dx, double beta)
{
	for (int i = 0; i < M; i++)
		F[i] = F1[i](x[0] + beta * dx[0], x[1] + beta * dx[1]);
}

void CalculateMatrix3(double** A, double* x)
{
	for (int i = 0; i < M; i++)
		for (int j = 0; j < N; j++)
			A[i][j] = Afuncs3[i][j](x[0], x[1], x[2]);
}

void CalculateF3(double* F, double* x)
{
	for (int i = 0; i < M; i++)
		F[i] = -Ffuncs3[i](x[0], x[1], x[2]);
}

void CalculateF3(double* F, double* x, double* dx, double beta)
{
	for (int i = 0; i < M; i++)
		F[i] = Ffuncs3[i](x[0] + beta * dx[0], x[1] + beta * dx[1], x[2] + beta * dx[2]);
}

void ProcessSLAE1(double** A, double* F)
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

void ProcessSLAE2(double** A, double* maxCol, int* indices, double* F)
{
	for (int i = 0; i < N; i++)
		indices[i] = i;

	for (int i = 0; i < N; i++)
	{
		double max = abs(A[0][i]);
		for (int k = 1; k < M; k++)
		{
			if (max < abs(A[k][i]))
				max = abs(A[k][i]);
		}

		maxCol[i] = max;
	}

	for (int i = 0; i < N; i++)
		for (int j = i + 1; j < N; j++)
		{
			if (maxCol[i] < maxCol[j])
			{
				swap(maxCol[i], maxCol[j]);
				swap(indices[i], indices[j]);
				for (int k = 0; k < M; k++)
					swap(A[k][i], A[k][j]);
			}
		}
}

void Multiply(double** A, double* vec, double* res)
{
	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			res[i] = A[i][j] * vec[j];
}

void Solve(int size, double** A, double* x, double* F)
{
	for (int k = 0; k < size - 1; k++)
	{
		for (int i = k + 1; i < size; i++)
		{
			double t = A[i][k] / A[k][k];
			F[i] -= t * F[k];
			for (int j = k + 1; j < size; j++)
				A[i][j] -= t * A[k][j];
		}
	}

	x[size - 1] = F[size - 1] / A[size - 1][size - 1];
	for (int k = size - 2; k >= 0; k--)
	{
		double sum = 0;
		for (int j = k + 1; j < size; j++)
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