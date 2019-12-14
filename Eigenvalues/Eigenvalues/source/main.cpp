#include "Matrix.h"

int main()
{
	Matrix A = { };
	int maxiter;
	double eps;
	ReadMatrix(A, maxiter, eps);
	double* x0 = new double[A.N];
	for (int i = 0; i < A.N; i++)
		x0[i] = 1;

	double* x1 = new double[A.N];
	double maxLambda = MaxEigenValue(A, x0, x1, maxiter, eps);
	double minLambda = MinEigenValue(A, x0, x1, maxiter, eps);
	return 0;
}