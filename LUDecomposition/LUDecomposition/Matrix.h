#pragma once

#include <fstream>
#include <iostream>

typedef double real;
typedef double realScal;

void InputSize(int&, int&);													// Entering A, AL, AU dimensions
void Input(int&, int&, real *, real *, real *, real *, int *);					// Entering matrix A in profile format
void BuildLU(int&, real *, real *, real *, int *);								// LU-decomposition
void Compute(int&, real *, real *, real *, int *, real *);						// Solution of a system of linear equations
void Output(int&, int&, real *, real *, real *, real *);						// The output of the matrices L and U in the profile format, the output of the vector x
void HilbertMatrix(int&, real *, real *, real *, int *);						// Create a Hilbert matrix of size N
void ToTightFormat(int&, real *, real *, real *, int *, real **);				// Convert matrix from profile format to tight format
void Gauss(int&, real*, real*, real**);										// Gauss method
void Multiply(int&, real*, real*, real*, int*, real *, real *);				// Multiplication matrix by vector