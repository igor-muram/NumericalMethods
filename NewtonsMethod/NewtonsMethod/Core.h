#pragma once
#include <functional>
#include <vector>

using namespace std;

const int N = 2;
const int M = 2;

void CalculateMatrix(double** A, double* x);
void CalculateF(double* F, double* x);
void CalculateF(double* F, double* x, double* dx, double beta);

void CalculateMatrix3(double** A, double* x);
void CalculateF3(double* F, double* x);
void CalculateF3(double* F, double* x, double* dx, double beta);

void ProcessSLAE1(double** A, double* F);
void ProcessSLAE2(double** A, double* maxCol, int* indices, double* F);
void Multiply(double** A, double* vec, double* res);
void Solve(int size, double** A, double* x, double* F);
double Norm(double* vec);