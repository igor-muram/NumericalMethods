#pragma once
#include <functional>
#include <vector>

using namespace std;

const int N = 2;

void CalculateMatrix(double** A, double* x);
void CalculateF(double* F, double* x);
void CalculateF(double* F, double* x, double* dx, double beta);
void Multiply(double** A, double* vec, double* res);
void Solve(double** A, double* x, double* F);
double Norm(double* vec);