#pragma once

#include <vector>

#include "Matrix.h"
#include "IntervalBuilder.h"
#include "FEMInfo.h"

void SetBoundary(ProfileMatrix& A, std::vector<double>& b, IntervalBuilder X, IntervalBuilder Y, IntervalBuilder Z)
{
	int xCount = X.Count(), yCount = Y.Count(), zCount = Z.Count();

	for (int i = 0; i < yCount; i++)
	{
		for (int j = 0; j < xCount; j++)
		{
			int line = i * xCount + j;

			A(2 * line, 2 * line) = 1.0e+50;
			b[2 * line] = 1.0e+50 * usin(X[j], Y[i], Z[0]);

			A(2 * line + 1, 2 * line + 1) = 1.0e+50;
			b[2 * line + 1] = 1.0e+50 * ucos(X[j], Y[i], Z[0]);
		}
	}

	for (int i = 0; i < zCount; i++)
	{
		for (int j = 0; j < xCount; j++)
		{
			int line = i * yCount * xCount + j;

			A(2 * line, 2 * line) = 1.0e+50;
			b[2 * line] = 1.0e+50 * usin(X[j], Y[0], Z[i]);

			A(2 * line + 1, 2 * line + 1) = 1.0e+50;
			b[2 * line + 1] = 1.0e+50 * ucos(X[j], Y[0], Z[i]);
		}
	}

	for (int i = 0; i < zCount; i++)
	{
		for (int j = 0; j < yCount; j++)
		{
			int line = i * yCount * xCount + j * xCount + xCount - 1;

			A(2 * line, 2 * line) = 1.0e+50;
			b[2 * line] = 1.0e+50 * usin(X.Back(), Y[j], Z[i]);

			A(2 * line + 1, 2 * line + 1) = 1.0e+50;
			b[2 * line + 1] = 1.0e+50 * ucos(X.Back(), Y[j], Z[i]);
		}
	}

	for (int i = 0; i < zCount; i++)
	{
		for (int j = 0; j < xCount; j++)
		{
			int line = i * yCount * xCount + (yCount - 1) * xCount + j;

			A(2 * line, 2 * line) = 1.0e+50;
			b[2 * line] = 1.0e+50 * usin(X[j], Y.Back(), Z[i]);

			A(2 * line + 1, 2 * line + 1) = 1.0e+50;
			b[2 * line + 1] = 1.0e+50 * ucos(X[j], Y.Back(), Z[i]);
		}
	}

	for (int i = 0; i < zCount; i++)
	{
		for (int j = 0; j < yCount; j++)
		{
			int line = i * yCount * xCount + j * xCount;

			A(2 * line, 2 * line) = 1.0e+50;
			b[2 * line] = 1.0e+50 * usin(X[0], Y[j], Z[i]);

			A(2 * line + 1, 2 * line + 1) = 1.0e+50;
			b[2 * line + 1] = 1.0e+50 * ucos(X[0], Y[j], Z[i]);
		}
	}

	for (int i = 0; i < yCount; i++)
	{
		for (int j = 0; j < xCount; j++)
		{
			int line = (zCount - 1) * xCount * yCount + i * xCount + j;

			A(2 * line, 2 * line) = 1.0e+50;
			b[2 * line] = 1.0e+50 * usin(X[j], Y[i], Z.Back());

			A(2 * line + 1, 2 * line + 1) = 1.0e+50;
			b[2 * line + 1] = 1.0e+50 * ucos(X[j], Y[i], Z.Back());
		}
	}
}

void SetBoundary(SparseMatrix& A, std::vector<double>& b, IntervalBuilder X, IntervalBuilder Y, IntervalBuilder Z)
{
	int xCount = X.Count(), yCount = Y.Count(), zCount = Z.Count();

	for (int i = 0; i < yCount; i++)
	{
		for (int j = 0; j < xCount; j++)
		{
			int line = i * xCount + j;

			A(2 * line, 2 * line) = 1.0e+50;
			b[2 * line] = 1.0e+50 * usin(X[j], Y[i], Z[0]);

			A(2 * line + 1, 2 * line + 1) = 1.0e+50;
			b[2 * line + 1] = 1.0e+50 * ucos(X[j], Y[i], Z[0]);
		}
	}

	for (int i = 0; i < zCount; i++)
	{
		for (int j = 0; j < xCount; j++)
		{
			int line = i * yCount * xCount + j;

			A(2 * line, 2 * line) = 1.0e+50;
			b[2 * line] = 1.0e+50 * usin(X[j], Y[0], Z[i]);

			A(2 * line + 1, 2 * line + 1) = 1.0e+50;
			b[2 * line + 1] = 1.0e+50 * ucos(X[j], Y[0], Z[i]);
		}
	}

	for (int i = 0; i < zCount; i++)
	{
		for (int j = 0; j < yCount; j++)
		{
			int line = i * yCount * xCount + j * xCount + xCount - 1;

			A(2 * line, 2 * line) = 1.0e+50;
			b[2 * line] = 1.0e+50 * usin(X.Back(), Y[j], Z[i]);

			A(2 * line + 1, 2 * line + 1) = 1.0e+50;
			b[2 * line + 1] = 1.0e+50 * ucos(X.Back(), Y[j], Z[i]);
		}
	}

	for (int i = 0; i < zCount; i++)
	{
		for (int j = 0; j < xCount; j++)
		{
			int line = i * yCount * xCount + (yCount - 1) * xCount + j;

			A(2 * line, 2 * line) = 1.0e+50;
			b[2 * line] = 1.0e+50 * usin(X[j], Y.Back(), Z[i]);

			A(2 * line + 1, 2 * line + 1) = 1.0e+50;
			b[2 * line + 1] = 1.0e+50 * ucos(X[j], Y.Back(), Z[i]);
		}
	}

	for (int i = 0; i < zCount; i++)
	{
		for (int j = 0; j < yCount; j++)
		{
			int line = i * yCount * xCount + j * xCount;

			A(2 * line, 2 * line) = 1.0e+50;
			b[2 * line] = 1.0e+50 * usin(X[0], Y[j], Z[i]);

			A(2 * line + 1, 2 * line + 1) = 1.0e+50;
			b[2 * line + 1] = 1.0e+50 * ucos(X[0], Y[j], Z[i]);
		}
	}

	for (int i = 0; i < yCount; i++)
	{
		for (int j = 0; j < xCount; j++)
		{
			int line = (zCount - 1) * xCount * yCount + i * xCount + j;

			A(2 * line, 2 * line) = 1.0e+50;
			b[2 * line] = 1.0e+50 * usin(X[j], Y[i], Z.Back());

			A(2 * line + 1, 2 * line + 1) = 1.0e+50;
			b[2 * line + 1] = 1.0e+50 * ucos(X[j], Y[i], Z.Back());
		}
	}
}
