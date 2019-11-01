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
	double* r = new double[A.N];
	double* z = new double[A.N];
	double* Ax = new double[A.N];
	double* p = new double[A.N];

	double* x = new double[A.N];
	for (int i = 0; i < A.N; i++)
		x[i] = 1;

	Multiply(A, x, f);
	for (int i = 0; i < A.N; i++)
		x[i] = 0;

	int iter = LOS(A, x, f, r, z, p, Ax, maxiter, eps);

	std::cout << "Number of iterations (LOS): " << iter << std::endl;
	for (int i = 0; i < A.N; i++)
		std::cout << x[i] << " ";
	std::cout << std::endl << std::endl;

	iter = CGM(A, x, f, r, z, Ax, maxiter, eps);

	std::cout << "Number of iterations (CGM): " << iter << std::endl;
	for (int i = 0; i < A.N; i++)
		std::cout << x[i] << " ";
	std::cout << std::endl << std::endl;

	system("pause");
	return 0;
}