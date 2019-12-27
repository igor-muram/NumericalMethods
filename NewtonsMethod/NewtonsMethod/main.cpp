#include "Core.h"
#include <fstream>
#include <vector>
#include <string>

using namespace std;

struct Point
{
	Point(double x, double y) : x(x), y(y) {  }
	double x, y;
};

struct AuxVectors
{
	double* maxCol = nullptr;
	double* dx = nullptr;
	int* indices = nullptr;
};

vector<Point> points;

void WritePoints(string filename)
{
	ofstream out(filename);
	out << points.size() << endl;
	for (auto point : points)
		out << point.x << "\t" << point.y << endl;

	out.close();
}

void ReadInitialValues(double* x0, double& eps, int& maxiter, double& delta)
{
	ifstream in("values.txt");
	for (int i = 0; i < N; i++)
		in >> x0[i];

	in >> eps >> maxiter >> delta;

	in.close();
}

void Newton1(double** A, double* F, double* x0, AuxVectors& aux, double eps, int maxiter, double delta)
{
	double* dx = aux.dx;
	double beta = 1;
	CalculateF(F, x0);

	points.emplace_back(x0[0], x0[1]);
	printf_s("x: %f\ty: %f\tbeta: %f\n", x0[0], x0[1], 1.0);

	for (int i = 0; i < maxiter && Norm(F) >= eps && beta >= delta; i++)
	{
		beta = 1;
		CalculateF(F, x0);
		CalculateMatrix(A, x0);
		Solve(N, A, dx, F);

		double norm = Norm(F);

		CalculateF(F, x0, dx, beta);
		while (Norm(F) >= norm && beta >= delta)
		{
			beta /= 2;
			CalculateF(F, x0, dx, beta);
		}

		if (beta >= delta)
		{
			for (int i = 0; i < N; i++)
				x0[i] += beta * dx[i];
		}

		points.emplace_back(x0[0], x0[1]);
		printf_s("x: %f\ty: %f\tbeta: %f\n", x0[0], x0[1], beta);
	}
}

void Newton2(double** A, double* F, double* x0, AuxVectors& aux, double eps, int maxiter, double delta)
{
	double* dx = aux.dx;
	double beta = 1;
	CalculateF(F, x0);

	points.emplace_back(x0[0], x0[1]);
	printf_s("x: %f\ty: %f\tbeta: %f\n", x0[0], x0[1], 1.0);

	for (int i = 0; i < maxiter && Norm(F) >= eps && beta >= delta; i++)
	{
		beta = 1;
		CalculateF(F, x0);
		CalculateMatrix(A, x0);
		ProcessSLAE1(A, F);
		Solve(N, A, dx, F);

		double norm = Norm(F);

		CalculateF(F, x0, dx, beta);
		while (Norm(F) >= norm && beta >= delta)
		{
			beta /= 2;
			CalculateF(F, x0, dx, beta);
		}

		for (int i = 0; i < N; i++)
			x0[i] += beta * dx[i];

		points.emplace_back(x0[0], x0[1]);
		printf_s("x: %f\ty: %f\tbeta: %f\n", x0[0], x0[1], beta);
	}
}

void Newton3(double** A, double* F, double* x0, AuxVectors& aux, double eps, int maxiter, double delta)
{
	double* maxCol = aux.maxCol;
	int* indices = aux.indices;
	double* dx = aux.dx;
	double beta = 1;
	CalculateF3(F, x0);

	printf_s("x: %f\ty: %f\tbeta: %f\n", x0[0], x0[1], x0[2], 1.0);

	for (int i = 0; i < maxiter && Norm(F) >= eps && beta >= delta; i++)
	{
		beta = 1;
		CalculateF3(F, x0);
		CalculateMatrix3(A, x0);
		ProcessSLAE2(A, maxCol, indices, F);
		Solve(M, A, dx, F);

		for (int i = M; i < N; i++)
			dx[i] = 0;

		for (int i = 0; i < M; i++)
			swap(dx[i], dx[indices[i]]);

		double norm = Norm(F);
		CalculateF3(F, x0, dx, beta);
		while (Norm(F) >= norm && beta >= delta)
		{
			beta /= 2;
			CalculateF3(F, x0, dx, beta);
		}

		if (beta >= delta)
		{
			for (int i = 0; i < N; i++)
				x0[i] += beta * dx[i];
		}

		printf_s("x: %f\ty: %f\tbeta: %f\n", x0[0], x0[1], x0[2], beta);
	}
}

void main()
{
	double* x0 = new double[N];
	double** A = new double* [M];
	for (int i = 0; i < M; i++)
		A[i] = new double[N];

	double* F = new double[M];
	double eps = 0, delta = 0;
	int maxiter = 0;
	ReadInitialValues(x0, eps, maxiter, delta);

	AuxVectors aux;
	aux.dx = new double[N];

	if (M == N)
	{
		Newton1(A, F, x0, aux, eps, maxiter, delta);
	}
	else if (M > N)
	{
		Newton2(A, F, x0, aux, eps, maxiter, delta);
	}
	else
	{
		aux.indices = new int[N];
		aux.maxCol = new double[N];
		Newton3(A, F, x0, aux, eps, maxiter, delta);
	}

	WritePoints("C:/input/points.txt");
	system("pause");
}