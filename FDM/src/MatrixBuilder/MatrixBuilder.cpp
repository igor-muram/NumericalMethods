#include "MatrixBuilder.h"

void BuildMatrix(
	SLAE::Matrix& A,
	vector<double>& b,
	vector<vector<int>>& areas,
	vector<Interval>& intervalsX,
	vector<Interval>& intervalsY,
	vector<BoundaryCondition>& conds,
	vector<double>& x,
	vector<double>& y,
	int kx, int ky)
{
	for (int iy = 1; iy < ky - 1; iy++)
	{
		int areaJ = IntervalNo(intervalsY, iy);

		for (int ix = 1; ix < kx - 1; ix++)
		{
			int areaI = IntervalNo(intervalsX, ix);

			if (areas[areaJ][areaI] != 0 && !IsOnBorder(conds, ix, iy))
			{
				int row = iy * kx + ix;

				double hx = x[ix + 1] - x[ix];
				double hxPrev = x[ix] - x[ix - 1];

				double hy = y[iy + 1] - y[iy];
				double hyPrev = y[iy] - y[iy - 1];

				A.D[row] = 2.0 / (hx * hxPrev) + 2.0 / (hy * hyPrev) + gamma;
				A.L1[row - 1] = -2.0 / (hxPrev * (hxPrev + hx));
				A.U1[row] = -2.0 / (hx * (hxPrev + hx));
				A.L2[row - kx] = -2.0 / (hyPrev * (hyPrev + hy));
				A.U2[row] = -2.0 / (hy * (hyPrev + hy));

				b[row] = f(x[ix], y[iy]);
			}
		}
	}

}

void BoundaryConditions(
	SLAE::Matrix& A,
	vector<double>& b,
	vector<double>& x,
	vector<double>& y,
	vector<BoundaryCondition>& conds)
{
	int kx = x.size();

	for (auto c : conds)
	{
		if (c.xBegin == c.xEnd)
		{
			int ix = c.xBegin;
			int begin = c.yBegin;
			int end = c.yEnd;

			for (int iy = begin; iy <= end; iy++)
			{
				int row = iy * kx + ix;
				A.D[row] = 1.0;
				b[row] = borderFuncs[c.functionNo](x[ix], y[iy]);
			}
		}
		else
		{
			int iy = c.yBegin;
			int begin = c.xBegin;
			int end = c.xEnd;

			for (int ix = begin; ix <= end; ix++)
			{
				int row = iy * kx + ix;
				A.D[row] = 1.0;
				b[row] = borderFuncs[c.functionNo](x[ix], y[iy]);
			}
		}
	}
}

int IntervalNo(vector<Interval>& intervals, int index)
{
	int areaI = 0;
	bool isFound = false;
	for (int i = 0; i < intervals.size() && !isFound; i++)
		if (index <= intervals[i].endN && index >= intervals[i].beginN)
		{
			areaI = i;
			isFound = true;
		}

	return areaI;
}

bool IsOnBorder(vector<BoundaryCondition>& conds, int ix, int iy)
{
	bool isOnBorder = false;
	int size = conds.size();

	for (int i = 0; i < size && !isOnBorder; i++)
	{
		if (ix >= conds[i].xBegin && ix <= conds[i].xEnd && iy >= conds[i].yBegin && iy <= conds[i].yEnd)
			isOnBorder = true;
	}

	return isOnBorder;
}