#include "MeshBuilder/MeshBuilder.h"
#include "MatrixBuilder/MatrixBuilder.h"

#include <SLAE/SLAE.h>

#include <iostream>
#include <vector>

using namespace std;

struct point {
	point(double x, double y) : x(x), y(y) {}
	double x, y;

	bool operator==(point a) {
		if (abs(a.x - x) < 1.0e-6 && abs(a.y - y) < 1.0e-6)
			return true;
		else
			return false;
	}
};

vector<point> points = {
		point(0.125, 0.125),
		point(0.125, 0.25),
		point(0.125, 0.375),
		point(0.125, 0.5),
		point(0.125, 0.625),
		point(0.125, 0.75),
		point(0.125, 0.875),

		point(0.25, 0.125),
		point(0.25, 0.25),
		point(0.25, 0.375),
		point(0.25, 0.5),
		point(0.25, 0.625),
		point(0.25, 0.75),
		point(0.25, 0.875),

		point(0.375, 0.125),
		point(0.375, 0.25),
		point(0.375, 0.375),
		point(0.375, 0.5),
		point(0.375, 0.625),
		point(0.375, 0.75),
		point(0.375, 0.875),

		point(0.5, 0.125),
		point(0.5, 0.25),
		point(0.5, 0.375),
		point(0.5, 0.5),
		point(0.5, 0.625),
		point(0.5, 0.75),
		point(0.5, 0.875),

		point(0.625, 0.125),
		point(0.625, 0.25),
		point(0.625, 0.375),
		point(0.625, 0.5),
		point(0.625, 0.625),
		point(0.625, 0.75),
		point(0.625, 0.875),

		point(0.75, 0.125),
		point(0.75, 0.25),
		point(0.75, 0.375),
		point(0.75, 0.5),
		point(0.75, 0.625),
		point(0.75, 0.75),
		point(0.75, 0.875),

		point(0.875, 0.125),
		point(0.875, 0.25),
		point(0.875, 0.375),
		point(0.875, 0.5),
		point(0.875, 0.625),
		point(0.875, 0.75),
		point(0.875, 0.875)
};


void PrintPoints(string filename, vector<double>& x, vector<double>& ix, vector<double>& iy, int kx, int ky)
{
	ofstream out(filename);

	for (int i = 0; i < kx; i++)
	{
		for (int j = 0; j < ky; j++)
		{
			bool flag = false;

			point p = point(ix[i], iy[j]);
			for (int k = 0; k < 49 && !flag; k++)
			{
				if (points[k] == p)
					flag = true;
			}

			if (flag)
				out << x[j * kx + i] << ';' << endl;
		}
	}

	out.close();
}

void PrintTest(string filename, vector<double>& x, int kx)
{
	ofstream out(filename);

	int size = x.size();
	int rowCount = size / kx;
	for (int i = 0; i < rowCount; i++)
	{
		for (int j = 0; j < kx; j++)
			out << x[i * kx + j] << ';';

		out << endl;
	}

	out.close();
}

int main()
{
	vector<Interval> intervalsX, intervalsY;
	ReadIntervals("input/intervalsX.txt", intervalsX);
	ReadIntervals("input/intervalsY.txt", intervalsY);

	vector<vector<int>> areas;
	ReadAreaMatrix("input/areaMatrix.txt", areas);

	vector<BoundaryCondition> conds;
	ReadBoundaryConds("input/boundary.txt", conds);

	int kx = CountNodes(intervalsX);
	int ky = CountNodes(intervalsY);

	IntervalNumbering(intervalsX);
	IntervalNumbering(intervalsY);

	BoundaryCondsNumbering(intervalsX, intervalsY, conds);

	vector<double> x(kx), y(ky);
	BuildMesh(intervalsX, x);
	BuildMesh(intervalsY, y);

	SLAE::Matrix A;
	vector<double> b(kx * ky);

	SLAE::InitMatrix(A, kx * ky, kx);
	BuildMatrix(A, b, areas, intervalsX, intervalsY, conds, x, y, kx, ky);
	BoundaryConditions(A, b, x, y, conds);

	vector<double> x0(kx * ky), r(kx * ky);
	int k = SLAE::Zeidel(A, b, x0, r, 1.7);

	PrintTest("result.csv", x0, kx);
	return 0;
}