#pragma once

#include <vector>
#include <iostream>

#include "FEMInfo.h"

using namespace std;

class SLAEBuilder
{
public:
	SLAEBuilder(int nodeCount, FEIterator begin, FEIterator end) : nodeCount(nodeCount), begin(begin), end(end)
	{ }

	void BuildGlobal(Matrix& A, vector<double>& q, double dt)
	{
		for (FEIterator iter = begin; iter != end; iter++)
		{
			FiniteElement e = *iter;
			vector<vector<double>> local = buildLocal(e, q, dt);

			for (int i = 0; i < 3; i++)
				for (int j = 0; j < 3; j++)
					A(e.v[i], e.v[j]) += local[i][j];
		}
	}

	void BuildGlobalB(vector<double>& b, vector<double>& q, double t, double dt)
	{
		for (FEIterator iter = begin; iter != end; iter++)
		{
			FiniteElement e = *iter;
			vector<double> local = buildLocalB(e, q, t, dt);

			for (int i = 0; i < 3; i++)
				b[e.v[i]] += local[i];
		}
	}

	void Boundary(Matrix& A, vector<double>& b, double u1, double u2)
	{
		A(0, 0) = 1.0e+50;
		b[0] = u1 * 1.0e+50;

		A(nodeCount - 1, nodeCount - 1) = 1.0e+50;
		b[nodeCount - 1] = u2 * 1.0e+50;
	}

private:
	FEIterator begin, end;
	int nodeCount;

	vector<vector<double>> buildLocal(FiniteElement e, vector<double>& qGlobal, double dt)
	{
		vector<vector<double>> local(3, vector<double>(3, 0));
		vector<double> q;
		double h = e.end - e.begin;

		for (auto v : e.v)
			q.push_back(qGlobal[v]);

		double l1 = lambda(e.begin, e.begin, h, q);
		double l2 = lambda(e.begin, (e.begin + e.end) / 2, h, q);
		double l3 = lambda(e.begin, e.end, h, q);

		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				local[i][j] = (l1 * G1[i][j] + l2 * G2[i][j] + l3 * G3[i][j]) / h + h * sigma * M[i][j] / dt;

		return local;
	}

	vector<double> buildLocalB(FiniteElement e, vector<double>& qGlobal, double t, double dt)
	{
		vector<double> f, q;
		f.push_back(F(e.begin, t));
		f.push_back(F((e.begin + e.end) / 2, t));
		f.push_back(F(e.end, t));

		for (auto v : e.v)
			q.push_back(qGlobal[v]);

		vector<double> b(3, 0.0);
		double h = e.end - e.begin;

		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				b[i] += h * (f[j] * M[i][j] + sigma * q[j] * M[i][j] / dt);

		return b;
	}

	double lambda(double b, double x, double h, vector<double> q)
	{
		double du1 = basisDirs[0](x - b, h) * q[0];
		double du2 = basisDirs[1](x - b, h) * q[1];
		double du3 = basisDirs[2](x - b, h) * q[2];
		return Lamdba(x, du1 + du2 + du3);
	}
};