#pragma once
#include <fstream>
#include <iostream>
#include <string>

using namespace std;

typedef double real;

const int maxiter = 10000;
const real eps = 1e-12;

struct Matrix
{
	real *D, *L1, *L2, *U1, *U2;
	int N, m;
};

void ReadMatrix(Matrix&);										// Entering matrix
void ReadX0(real *, int&);										// Entering initial approximation
void ReadF(real *, int&);										// Entering F
void Output(int&, real *);										// The output of the vector x
real Multiply(Matrix&, real *, int);							// Multiplication matrix by vector
real Norm(real *, int&);										// Norm of vector
int Jacobi(Matrix&, real *, real *, real *, real *, real);		// Jacobi Method
int Zeidel(Matrix&, real *, real *, real *, real);				// Zeidel Method