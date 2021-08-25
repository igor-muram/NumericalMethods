#pragma once
#include "stdafx.h"
#include "FEMInfo.h"
#include "Matrix.h"

class SLAEBuilder
{
public:
	SLAEBuilder(FEIterator begin, FEIterator end) : begin(begin), end(end) {}

	void BuildSparseMatrix(SparseMatrix& A)
	{
		for (FEIterator iter = begin; iter < end; iter++)
		{
			FE e = *iter;
			std::vector<std::vector<double>> local = buildLocal(iter);
			
			for (int i = 0; i < 8; i++)
				for (int j = 0; j < 8; j++)
				{
					int Ai = e.verts[i];
					int Aj = e.verts[j];
					A(2 * Ai, 2 * Aj) += local[2 * i][2 * j];
					A(2 * Ai, 2 * Aj + 1) += local[2 * i][2 * j + 1];
					A(2 * Ai + 1, 2 * Aj) += local[2 * i + 1][2 * j];
					A(2 * Ai + 1, 2 * Aj + 1) += local[2 * i + 1][2 * j + 1];
				}
		}
	}

	void BuildProfileMatrix(ProfileMatrix& A)
	{
		for (FEIterator iter = begin; iter < end; iter++)
		{
			std::vector<std::vector<double>> local = buildLocal(iter);

			for (int i = 0; i < 8; i++)
				for (int j = 0; j < 8; j++)
				{
					int Ai = iter->verts[i];
					int Aj = iter->verts[j];
					A(2 * Ai, 2 * Aj) += local[2 * i][2 * j];
					A(2 * Ai, 2 * Aj + 1) += local[2 * i][2 * j + 1];
					A(2 * Ai + 1, 2 * Aj) += local[2 * i + 1][2 * j];
					A(2 * Ai + 1, 2 * Aj + 1) += local[2 * i + 1][2 * j + 1];
				}
		}
	}

	void BuildB(std::vector<double>& b)
	{
		for (FEIterator iter = begin; iter != end; iter++)
		{
			std::vector<double> local = buildLocalB(iter);

			for (int i = 0; i < 8; i++)
			{
				int bi = iter->verts[i];
				b[2 * bi] += local[2 * i];
				b[2 * bi + 1] += local[2 * i + 1];
			}
		}
	}

private:
	FEIterator begin, end;

private:
	std::vector<std::vector<double>> buildLocal(FEIterator e)
	{
		std::vector<std::vector<double>> local(16, std::vector<double>(16, 0));
		double hx = e->x[1] - e->x[0];
		double hy = e->y[1] - e->y[0];
		double hz = e->z[1] - e->z[0];

		for (int i = 0; i < 8; i++)
			for (int j = 0; j < 8; j++)
			{
				double g = G(i, j, hx, hy, hz);
				double m = M(i, j, hx, hy, hz);

				double p = lambda * g - w * w * chi * m;
				double c = w * sigma * m;
				local[2 * i][2 * j] = p;
				local[2 * i][2 * j + 1] = -c;
				local[2 * i + 1][2 * j] = c;
				local[2 * i + 1][2 * j + 1] = p;
			}

		return local;
	}

	std::vector<double> buildLocalB(FEIterator e)
	{
		double hx = e->x[1] - e->x[0];
		double hy = e->y[1] - e->y[0];
		double hz = e->z[1] - e->z[0];

		std::vector<double> fs, fc;
		for (int i = 0; i < 2; i++)
			for (int j = 0; j < 2; j++)
				for (int k = 0; k < 2; k++)
				{
					fs.push_back(fsin(e->x[k], e->y[j], e->z[i]));
					fc.push_back(fcos(e->x[k], e->y[j], e->z[i]));
				}

		std::vector<double> bs(8), bc(8);
		for (int i = 0; i < 8; i++)
			for (int j = 0; j < 8; j++)
			{
				bs[i] += fs[j] * M(i, j, hx, hy, hz);
				bc[i] += fc[j] * M(i, j, hx, hy, hz);
			}

		std::vector<double> b;

		for (int i = 0; i < 8; i++)
		{
			b.push_back(bs[i]);
			b.push_back(bc[i]);
		}

		return b;
	}
};