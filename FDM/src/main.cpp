#include "MeshBuilder/MeshBuilder.h"
#include "MatrixBuilder/MatrixBuilder.h"

#include <SLAE/SLAE.h>

#include <iostream>
#include <iomanip>
#include <vector>

using namespace std;

struct point {
	point(double x, double y) : x(x), y(y) {}
	double x, y;

	bool operator==(point a) {
		return (abs(a.x - x) < 1.0e-6 && abs(a.y - y) < 1.0e-6);
	}
};

vector<point> points = {
		point(0.24031880500416805, 0.24031880500416805),
		point(0.24031880500416805, 0.43257384900750251),
		point(0.24031880500416805, 0.58637788421017012),
		point(0.24031880500416805, 0.70942111237230421),
		point(0.24031880500416805, 0.80785569490201148),
		point(0.24031880500416805, 0.88660336092577730),
		point(0.24031880500416805, 0.94960149374478997),

		point(0.43257384900750251, 0.24031880500416805),
		point(0.43257384900750251, 0.43257384900750251),
		point(0.43257384900750251, 0.58637788421017012),
		point(0.43257384900750251, 0.70942111237230421),
		point(0.43257384900750251, 0.80785569490201148),
		point(0.43257384900750251, 0.88660336092577730),
		point(0.43257384900750251, 0.94960149374478997),

		point(0.58637788421017012, 0.24031880500416805),
		point(0.58637788421017012, 0.43257384900750251),
		point(0.58637788421017012, 0.58637788421017012),
		point(0.58637788421017012, 0.70942111237230421),
		point(0.58637788421017012, 0.80785569490201148),
		point(0.58637788421017012, 0.88660336092577730),
		point(0.58637788421017012, 0.94960149374478997),

		point(0.70942111237230421, 0.24031880500416805),
		point(0.70942111237230421, 0.43257384900750251),
		point(0.70942111237230421, 0.58637788421017012),
		point(0.70942111237230421, 0.70942111237230421),
		point(0.70942111237230421, 0.80785569490201148),
		point(0.70942111237230421, 0.88660336092577730),
		point(0.70942111237230421, 0.94960149374478997),

		point(0.80785569490201148, 0.24031880500416805),
		point(0.80785569490201148, 0.43257384900750251),
		point(0.80785569490201148, 0.58637788421017012),
		point(0.80785569490201148, 0.70942111237230421),
		point(0.80785569490201148, 0.80785569490201148),
		point(0.80785569490201148, 0.88660336092577730),
		point(0.80785569490201148, 0.94960149374478997),

		point(0.88660336092577730, 0.24031880500416805),
		point(0.88660336092577730, 0.43257384900750251),
		point(0.88660336092577730, 0.58637788421017012),
		point(0.88660336092577730, 0.70942111237230421),
		point(0.88660336092577730, 0.80785569490201148),
		point(0.88660336092577730, 0.88660336092577730),
		point(0.88660336092577730, 0.94960149374478997),

		point(0.94960149374478997, 0.24031880500416805),
		point(0.94960149374478997, 0.43257384900750251),
		point(0.94960149374478997, 0.58637788421017012),
		point(0.94960149374478997, 0.70942111237230421),
		point(0.94960149374478997, 0.80785569490201148),
		point(0.94960149374478997, 0.88660336092577730),
		point(0.94960149374478997, 0.94960149374478997)
};


void PrintPoints(string filename, vector<double>& x, vector<double>& ix, vector<double>& iy, int kx, int ky)
{
	ofstream out(filename);

	for (int i = 0; i < kx; i++)
		for (int j = 0; j < ky; j++)
		{
			bool flag = false;

			point p = point(ix[i], iy[j]);
			for (int k = 0; k < 49 && !flag; k++)
				if (points[k] == p)
					flag = true;

			if (flag)
				out << x[j * kx + i] << ';' << endl;
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

	/*ofstream out("x.txt");
	for (int i = 0; i < 49; i++)
	{

		out << setprecision(15) << points[i].x << endl;
	}

	out.close();*/

	return 0;

}