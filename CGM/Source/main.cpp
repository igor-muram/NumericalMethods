#include "Matrix/matrix.h"
#include "CGM/cgm.h"
#include "LOS/los.h"
#include <cstdio>
#include <ctime>

int main()
{
	Matrix A = { };
	Matrix LU = { };
	AuxVectors aux = { };
	int maxiter = 100000, choice = 11;
	double eps;
	ReadMatrix(A, maxiter, eps, choice);
	//A = HilbertMatrix(40);

	double* f = new double[A.N];
	double* x = new double[A.N];
	
	for (int i = 0; i < A.N; i++)
		x[i] = i + 1;

	Multiply(A, x, f);
	ReadX0(A.N, x);

	aux.Ax = new double[A.N];
	aux.r = new double[A.N];
	aux.z = new double[A.N];
	
	int iter = 0;
	double lastdiff = 0;
	double norm_f = sqrt(DotProduct(A.N, f, f));

	clock_t start_time, end_time;
	start_time = clock();
	switch (choice)
	{
	case 11:
		iter = CGM(A, x, f, aux, maxiter, eps, lastdiff);
		break;

	case 12:
		aux.LU = new double[A.N];
		aux.temp = new double[A.N];
		iter = CGM_diag1(A, x, f, aux, maxiter, eps, lastdiff);
		break;

	case 13:
		aux.LU = new double[A.N];
		aux.temp = new double[A.N];
		iter = CGM_diag2(A, x, f, aux, maxiter, eps, lastdiff);
		break;

	case 14:
		aux.LU = new double[A.N];
		aux.temp = new double[A.N];
		iter = CGM_diag3(A, x, f, aux, maxiter, eps, lastdiff);
		break;

	case 15:
		aux.temp = new double[A.N];
		LUFactorization(A, LU);
		iter = CGM_LU(A, x, f, LU, aux, maxiter, eps, lastdiff);
		break;

	case 21:
		aux.p = new double[A.N];
		iter = LOS(A, x, f, aux, maxiter, eps, lastdiff);
		break;

	case 22:
		aux.p = new double[A.N];
		aux.LU = new double[A.N];
		aux.temp = new double[A.N];
		iter = LOS_diag1(A, x, f, aux, maxiter, eps, lastdiff);
		break;

	case 23:
		aux.p = new double[A.N];
		aux.LU = new double[A.N];
		iter = LOS_diag2(A, x, f, aux, maxiter, eps, lastdiff);
		break;

	case 24:
		aux.p = new double[A.N];
		aux.LU = new double[A.N];
		aux.temp = new double[A.N];
		iter = LOS_diag3(A, x, f, aux, maxiter, eps, lastdiff);
		break;

	case 25:
		aux.p = new double[A.N];
		aux.temp = new double[A.N];
		LUFactorization(A, LU);
		iter = LOS_LU(A, x, f, LU, aux, maxiter, eps, lastdiff);
		break;

	default:
		printf_s("ERROR! Method is not specified.");
	}
	end_time = clock();

	printf_s("time: %.3f\n", (double)(end_time - start_time) / CLOCKS_PER_SEC);
	printf_s("iteration count: %d\n", iter);
	printf_s("iterative difference: %e\n", lastdiff);
	for (int i = 0; i < A.N; i++)
		printf_s("%.15f\n", x[i]);

	printf_s("\n");

	// Real difference
	double* r = new double[A.N];
	for (int i = 0; i < A.N; i++)
		r[i] = 0;

	Multiply(A, x, r);
	double nev = 0, dif = 0;
	for (int i = 0; i < A.N; i++)
	{
		dif = f[i] - r[i];
		nev += dif * dif;
	}
	nev = sqrt(nev) / norm_f;
	printf_s("real difference: %e\n", nev);

	system("pause");
	return 0;
}