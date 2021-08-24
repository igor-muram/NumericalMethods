#include <iostream>

#include "BCG.h"
#include "Matrix.h"

int main()
{
	Matrix A;
	A.N = 4;
	A.DI = { 1, 8, 7, 2 };
	A.IA = { 0, 0, 1, 3, 6 };
	A.JA = { 0, 0, 1, 0, 1, 2 };
	A.AL = { 5, 3, 4, -2, 1, 8 };
	A.AU = { 1, 8, 2, 3, 5, 4 };

	vector<double> x = { 1.0, 1.0, 1.0, 1.0 };
	vector<double> b(4);
	vector<double> s = { 0.0, 0.0, 0.0, 0.0 };

	A.Multiply(x, b);

	BCG(A, s, b);

	std::cout << "Hello, World!" << std::endl;
	return 0;
}