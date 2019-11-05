#pragma once
#include <fstream>

using namespace std;

struct Matrix
{
	int N;
	double *DI, *AL, *AU;
	int *IA, *JA;
};

void ReadMatrix(Matrix& A, int& Aaxiter, double& eps, int& choice);
void ReadX0(int N, double *x0);
void ReadB(int N, double *b);
void Multiply(Matrix& A, double *vec, double *res);
void MultiplyT(Matrix& A, double *vec, double *res);
void LUFactorization(Matrix& A, Matrix& LU);
double DotProduct(int N, double *a, double *b);