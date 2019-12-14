#pragma once
#include <fstream>
#include <iostream>
#include <utility>


struct Matrix
{
	double* DI = nullptr;
	double* AL = nullptr;
	double* AU = nullptr;
	int* IA = nullptr;
	int N = 0;
};

void ReadMatrix(Matrix& A, int& maxiter, double& eps);			
void LUDecomposition(Matrix& A);
void Solve(Matrix& A, double* b, double* x);
void Multiply(Matrix& A, double* vec, double* res);
double Norm(int N, double* vector);
double MaxEigenValue(Matrix& A, double* x0, double* x1, int maxiter, double eps);
double MinEigenValue(Matrix& A, double* x0, double* x1, int maxiter, double eps);