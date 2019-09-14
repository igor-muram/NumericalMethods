#pragma once

#include <fstream>
#include <iostream>

typedef float real;

void InputSize(int&, int&);														// ���� ������� ������ A, AL, AU

void Input(int&, int&, real *, real *, real *, real *, int *);					// ���� ������� A � ���������� �������

void BuildLU(int&, real *, real *, real *, int *);								// LU-����������

void Compute(int&, real *, real *, real *, int *, real *);						// ������� ����

void Print(int&, int&, real *, real *, real *, real *);							// ����� ������ L � U � ���������� �������, ����� ������� x