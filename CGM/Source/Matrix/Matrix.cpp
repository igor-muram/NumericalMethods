#include "Matrix.h"					

void ReadMatrix(Matrix& A, int& Aaxiter, double& eps)
{
	std::ifstream in("input/kuslau.txt");
	in >> A.N >> Aaxiter >> eps;
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

void ReadB(int N, double* b)
{
	std::ifstream in("input/pr.txt");
	for (int i = 0; i < N; i++)
		in >> b[i];
	in.close();
}

void ReadX0(int N, double* x0)
{
	std::ifstream in("x0.txt");

	int i = 0;
	for (; in >> x0[i] && i < N; i++) {}

	for (; i < N; i++)
		x0[i] = 0;
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
		res[i] = vec[i] * di[i];

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

double DotProduct(int N, double* a, double* b)
{
	double sum = 0;
	for (int i = 0; i < N; i++)
		sum += a[i] * b[i];

	return sum;
}