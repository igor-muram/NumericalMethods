#include "MatrixBuilder.h"

void BuildMatrix(
	Matrix& A,
	vector<vector<int>>& areas,
	vector<Interval>& intervalsX,
	vector<Interval>& intervalsY,
	vector<double>& x,
	vector<double>& y,
	vector<double>& hx,
	vector<double>& hy,
	int kx, int ky)
{
	for (int iy = 1; iy < ky - 1; iy++)
	{
		int areaJ = IntervalNo(intervalsY, iy);

		for (int ix = 1; ix < kx - 1; ix++)
		{
			int areaI = IntervalNo(intervalsX, ix);

			if (areas[areaI][areaJ] != 0)
			{
				int row = iy * kx + ix;

				A.D[row] = 2.0 / (hx[ix] * hx[ix - 1]) + 2.0 / (hy[iy] * hy[iy - 1]) + gamma;
				A.L1[row - 1] = -2.0 / (hx[ix - 1] * (hx[ix - 1] + hx[ix]));
				A.U1[row + 1] = -2.0 / (hx[ix] * (hx[ix - 1] + hx[ix]));
				A.L2[row - (kx + 1)] = -2.0 / (hy[iy - 1] * (hy[iy - 1] + hy[iy]));
				A.U2[row + (kx + 1)] = -2.0 / (hy[iy] * (hy[iy - 1] + hy[iy]));
			}
		}
	}

}

int IntervalNo(vector<Interval>& intervals, int index)
{
	int areaI = 0;
	for (int i = 0; i < intervals.size(); i++)
		if (index <= intervals[i].endN && index >= intervals[i].beginN)
			areaI = i;

	return areaI;
}