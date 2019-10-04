#pragma once
#include <fstream>
#include <iostream>
#include <string>

using namespace std;

typedef float real;

void InputSize(int&, int&, int&);
void Input(int&, int&, real **, real *);
void Multiply(int&, int&, real **, int&, real *, real *);
void Output(int&, real *);