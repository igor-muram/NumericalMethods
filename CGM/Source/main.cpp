#include "Matrix/Matrix.h"
#include "CGM/cgm.h"
#include "LOS/los.h"
#include <iostream>

int main()
{
	Matrix A = { };
	int maxiter;
	double eps;
	ReadMatrix(A, maxiter, eps);

	double* f = new double[A.N];

	double* x0 = new double[A.N];
	for (int i = 0; i < A.N; i++)
		x0[i] = 0;

	double* r = new double[A.N];
	double* z = new double[A.N];

	double* x = new double[A.N];
	for (int i = 0; i < A.N; i++)
		x[i] = 1;

	Multiply(A, x, f);

	double* p = new double[A.N];
	int iter = LOS(A, x0, f, r, z, p, maxiter, eps, x);

	std::cout << "Number of iterations (LOS): " << iter << std::endl;
	for (int i = 0; i < A.N; i++)
		std::cout << x[i] << " ";
	std::cout << std::endl << std::endl;

	iter = CGM(A, x0, f, r, z, maxiter, eps, x);

	for (int i = 0; i < A.N; i++)
		x0[i] = 0;

	std::cout << "Number of iterations (CGM): " << iter << std::endl;
	for (int i = 0; i < A.N; i++)
		std::cout << x[i] << " ";
	std::cout << std::endl << std::endl;

	system("pause");
	return 0;
}