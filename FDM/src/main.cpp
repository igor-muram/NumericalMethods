#include "MeshBuilder/MeshBuilder.h"
#include "MatrixBuilder/Matrix.h"

#include <iostream>
#include <vector>

using namespace std;

int main()
{
	vector<Interval> intervalsX, intervalsY;
	ReadIntervals("input/intervalsX.txt", intervalsX);
	ReadIntervals("input/intervalsY.txt", intervalsY);

	vector<vector<int>> areas;
	ReadAreaMatrix("input/areaMatrix.txt", areas);

	int kx = CountBreakPoints(intervalsX);
	int ky = CountBreakPoints(intervalsY);

	IntervalNumbering(intervalsX, 1);
	IntervalNumbering(intervalsY, kx);

	vector<double> x(kx), y(ky);
	vector<double> hx(kx - 1), hy(ky - 1);

	BuildMesh(intervalsX, kx, x, hx);
	BuildMesh(intervalsY, ky, y, hy);

	Matrix A;
	InitMatrix(A, kx * ky, kx);

	return 0;
}