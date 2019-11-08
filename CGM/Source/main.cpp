#include "Matrix/matrix.h"
#include "CGM/cgm.h"
#include "LOS/los.h"
#include <iostream>

int main()
{
	Matrix A = { };
	AuxVectors aux = { };
	int maxiter, choice;
	double eps;
	ReadMatrix(A, maxiter, eps, choice);

	double* f = new double[A.N];
	double* x = new double[A.N];

	for (int i = 0; i < A.N; i++)
		x[i] = 1;

	Multiply(A, x, f);
	for (int i = 0; i < A.N; i++)
		x[i] = 0;

	aux.Ax = new double[A.N];
	aux.r = new double[A.N];
	aux.z = new double[A.N];

	switch (choice)
	{
	case 11:
		CGM(A, x, f, aux, maxiter, eps);
		break;

	case 12:
		aux.LU = new double[A.N];
		aux.temp = new double[A.N];
		CGM_diag1(A, x, f, aux, maxiter, eps);
		break;

	case 13:
		aux.LU = new double[A.N];
		aux.temp = new double[A.N];
		CGM_diag2(A, x, f, aux, maxiter, eps);
		break;

	case 14:
		aux.LU = new double[A.N];
		aux.temp = new double[A.N];
		CGM_diag3(A, x, f, aux, maxiter, eps);
		break;

	case 21:
		aux.p = new double[A.N];
		LOS(A, x, f, aux, maxiter, eps);
		break;

	case 22:
		aux.p = new double[A.N];
		aux.LU = new double[A.N];
		aux.temp = new double[A.N];
		LOS_diag1(A, x, f, aux, maxiter, eps);
		break;

	case 23:
		aux.p = new double[A.N];
		aux.LU = new double[A.N];
		LOS_diag2(A, x, f, aux, maxiter, eps);
		break;

	case 24:
		aux.p = new double[A.N];
		aux.LU = new double[A.N];
		aux.temp = new double[A.N];
		LOS_diag3(A, x, f, aux, maxiter, eps);
		break;

	default:
		std::cout << "ERROR! Method is not specified." << std::endl;
	}

	for (int i = 0; i < A.N; i++)
		std::cout << x[i] << " ";

	std::cout << std::endl << std::endl;

	system("pause");
	return 0;
}