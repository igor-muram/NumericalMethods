#include "Core.h"
#include <iostream>

using namespace std;

void ReadInitialValues(double* x0, double eps, int maxiter, double delta)
{
	ifstream in("values.txt");
	for (int i = 0; i < N; i++)
		in >> x0[i];

	in >> eps >> maxiter >> delta;
}

void Newton(double* x0, double eps, int maxiter, double delta)
{

}

void main()
{
	double* x0 = new double[N];
	double* A = new double[N * N];
	double* F = new double[N];
	double eps, delta;
	int maxiter;
	ReadInitialValues(x0, eps, maxiter, delta);
	Newton(x0, eps, maxiter, delta);
}