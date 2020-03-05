#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
using namespace std;

namespace SLAE 
{
	const int maxiter = 50000;
	const double eps = 1e-10;

	struct Matrix
	{
		vector<double> D, L1, L2, U1, U2;
		int N, m;
	};

	using TightMatrix = vector<vector<double>>;

	void InitMatrix(Matrix& A, int N, int m);
	void ToTight(Matrix& A, TightMatrix& TM);
	void PrintMatrix(Matrix& A);
	void PrintSLAEToCSV(string filename, Matrix& A, vector<double>& b);
	double Multiply(Matrix& M, vector<double>& vec, int i);
	double Norm(vector<double>& vec, int N);
	int Jacobi(Matrix& M, vector<double>& f, vector<double>& x0, vector<double>& x1, vector<double>& r, double w);
	int Zeidel(Matrix& M, vector<double>& f, vector<double>& x0, vector<double>& r, double w);
}
