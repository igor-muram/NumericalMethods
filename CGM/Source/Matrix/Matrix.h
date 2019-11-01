#pragma once
#include <fstream>

struct Matrix
{
	int N;
	double* DI;
	double* AL;
	double* AU;
	int* IA;
	int* JA;
};

void ReadMatrix(Matrix& A, int& Aaxiter, double& eps);
void ReadX0(int N, double* x0);
void ReadB(int N, double* b);
void Multiply(Matrix& A, double* vec, double* res);
double DotProduct(int N, double* a, double* b);