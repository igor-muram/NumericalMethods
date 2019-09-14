#pragma once

#include <fstream>
#include <iostream>

typedef float real;

void InputSize(int&, int&);														// Ввод размера матриц A, AL, AU

void Input(int&, int&, real *, real *, real *, real *, int *);					// Ввод матрицы A в профильном формате

void BuildLU(int&, real *, real *, real *, int *);								// LU-разложение

void Compute(int&, real *, real *, real *, int *, real *);						// Решение СЛАУ

void Print(int&, int&, real *, real *, real *, real *);							// Вывод матриц L и U в профильном формате, вывод вектора x