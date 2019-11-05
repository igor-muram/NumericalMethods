#include "Matrix/matrix.h"
#include "CGM/cgm.h"
#include "LOS/los.h"
#include "Memory/memory.h"
#include <iostream>

int main()
{
	Matrix A = { };
	int maxiter, choice;
	double eps;
	
	ReadMatrix(A, maxiter, eps, choice);

	Memory cache(A.N, 10);

	double *f;
	double *x;
	cache.Allocate(&f);
	cache.Allocate(&x);

	for (int i = 0; i < A.N; i++)
		x[i] = 1;

	Multiply(A, x, f);
	for (int i = 0; i < A.N; i++)
		x[i] = 0;

	int iter;
	if (choice == 1)
		iter = CGM_diag(A, x, f, cache, maxiter, eps);

	if (choice == 2)
		iter = LOS_diag(A, x, f, cache, maxiter, eps);

	system("pause");
	return 0;
}