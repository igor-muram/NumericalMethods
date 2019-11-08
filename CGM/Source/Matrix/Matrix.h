#pragma once
#include <fstream>

struct Matrix
{
	int N;
	double* DI = nullptr;
	double* AL = nullptr;
	double* AU = nullptr;
	int* IA = nullptr;
	int* JA = nullptr;
};

struct AuxVectors
{
	double* Ax = nullptr;
	double* r = nullptr;
	double* z = nullptr;
	double* p = nullptr;
	double* LU = nullptr;
	double* temp = nullptr;
};

void ReadMatrix(Matrix& A, int& Aaxiter, double& eps, int& choice);
void ReadX0(int N, double* x0);
void ReadB(int N, double* b);
void Multiply(Matrix& A, double* vec, double* res);
void MultiplyT(Matrix& A, double* vec, double* res);
void LUFactorization(Matrix& A, Matrix& LU);
double DotProduct(int N, double* a, double* b);