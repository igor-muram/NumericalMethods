#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>

#include "MeshBuilder.h"
#include "PortraitBuilder.h"
#include "SLAEBuilder.h"
#include "TimeMeshBuilder.h"
#include "CSV.h"
#include "Matrix.h"
#include "SLAESolver/SLAESolver.h"


std::vector<double> _ = { 0.0, 4.444444444444, 7.11111111111111111111, 10.666666666666666667, 14.22222222222222222, 16.0 };
std::vector<double> __ = { 0.0, 4.0156862745098039, 6.0235294117647058, 7.0274509803921568, 7.5294117647058822, 7.7803921568627450, 7.9058823529411768, 7.9686274509803923, 8.0 };

std::vector<double> _t = { 0.0, 0.5, 1.0, 1.5, 2.0 };
std::vector<double> __t = { 0.0, 1.0666666666666667, 1.6000000000000001, 1.8666666666666667, 2.0 };

double u(double x, double t)
{
	return x;
}

int main()
{
	// Build Mesh
	std::vector<Interval> intervals;
	std::ifstream in("input/intervals.txt");
	int count;
	in >> count;
	for (int i = 0; i < count; i++)
	{
		Interval interval;
		in >> interval;
		intervals.push_back(interval);
	}
	in.close();

	MeshBuilder mesh(intervals);

	// Build Time Mesh
	std::vector<Interval> time_intervals;
	in.open("input/time_intervals.txt");
	int time_count;
	in >> time_count;
	for (int i = 0; i < time_count; i++)
	{
		Interval interval;
		in >> interval;
		time_intervals.push_back(interval);
	}
	in.close();
	TimeMeshBuilder time(time_intervals);

	// Matrix
	Matrix A(mesh.GetNodeCount());

	// Build Matrix Portrait
	PortraitBuilder portrait(mesh.GetNodeCount(), mesh.Begin(), mesh.End());
	portrait.Build(A);
	A.AL.resize(A.IA.back());
	A.AU.resize(A.IA.back());

	SLAEBuilder slae(mesh.GetNodeCount(), mesh.Begin(), mesh.End());

	// Preparing for fixed-point iteration
	int n = mesh.GetNodeCount();
	vector<double> x, q0, q1(n, 0);

	// Find x coords
	FEIterator begin = mesh.Begin();
	FEIterator end = mesh.End();

	for (FEIterator i = begin; i != end; i++)
	{
		x.push_back(i->begin);
		x.push_back((i->begin + i->end) / 2);
	}
	x.push_back((end - 1)->end);

	// Find q0 
	TimeIterator time_begin = time.Begin();
	TimeIterator time_end = time.End();

	for (int i = 0; i < n; i++)
		q0.push_back(u(x[i], *time_begin));

	// CSV writer
	CSV csv(15, 100000);
	int col = 0;

	for (int i = 0, j = 0; i < n; i++)
		for (int k = 0; k < _.size(); k++)
			if (abs(x[i] - _[k]) < 1.0e-7)
			{
				csv(j + 1, col, x[i]);
				j++;
			}

	col += 2;

	for (int k = 0; k < _t.size(); k++)
	{
		if (abs(0.0 - _t[k]) < 1.0e-7)
		{
			csv(0, col, "t = ");
			csv(0, col, 0.0);
			for (int i = 0, j = 0; i < n; i++)
				for (int k = 0; k < _.size(); k++)
					if (abs(x[i] - _[k]) < 1.0e-7)
					{
						csv(j + 1, col, q0[i]);
						j++;
					}
			col += 2;
		}
	}

	// Start fixed-point iteration
	vector<double> b(n);
	vector<double> Aq(n);

	for (TimeIterator i = time_begin + 1; i != time_end; i++)
	{
		double t = *i;
		double tPrev = *(i - 1);
		double dt = t - tPrev;

		double eps = 10e-9;
		double diff = 1;

		double delta = 10e-9;
		double diff1 = 1;

		while (diff >= eps && diff1 >= delta)
		{
			slae.BuildGlobal(A, q0, dt);
			slae.BuildGlobalB(b, q0, t, dt);

			slae.Boundary(A, b, u(x[0], t), u(x[n - 1], t));

			LUDecomposition(A);
			Solve(A, q1, b);

			A.Clear();
			for (auto& bi : b)
				bi = 0;

			slae.BuildGlobal(A, q1, dt);
			slae.BuildGlobalB(b, q0, t, dt);

			Multiply(A, q1, Aq);

			Aq[0] = u(x[0], t);
			Aq[n - 1] = u(x[n - 1], t);

			b[0] = u(x[0], t);
			b[n - 1] = u(x[n - 1], t);

			diff = 0;
			double norm = 0;
			for (int i = 0; i < n; i++)
			{
				diff += (Aq[i] - b[i]) * (Aq[i] - b[i]);
				norm += b[i] * b[i];
			}
			diff = sqrt(diff / norm);

			diff1 = 0;
			double norm1 = 0;
			for (int i = 0; i < n; i++)
			{
				diff1 += (q1[i] - q0[i]) * (q1[i] - q0[i]);
				norm1 += q1[i] * q1[i];
			}
			diff1 = sqrt(diff1 / norm1);

			q0 = q1;

			for (int k = 0; k < _t.size(); k++)
			{
				if (abs(t - _t[k]) < 1.0e-7)
				{
					csv(0, col, "t = ");
					csv(0, col, t);
					int j = 0;
					for (int i = 0; i < n; i++)
						for (int k = 0; k < _.size(); k++)
							if (abs(x[i] - _[k]) < 1.0e-7)
							{
								csv(j + 1, col, q0[i]);
								j++;
							}

					csv(j + 1, col, diff);
					csv(j + 2, col, diff1);

					col += 2;
				}
			}

			A.Clear();
			for (auto& bi : b)
				bi = 0;
		}
	}

	csv.Write("result.csv");
	system("result.csv");

	return 0;
}