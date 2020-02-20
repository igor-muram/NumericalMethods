#include "MeshBuilder.h"

void ReadIntervals(string filename, vector<Interval>& intervals)
{
	ifstream in(filename);

	int count;
	in >> count;

	for (int i = 0; i < count; i++)
	{
		Interval interval;
		in >> interval.begin >> interval.end >> interval.n;
		intervals.push_back(interval);
	}

	in.close();
}

void ReadAreaMatrix(string filename, vector<vector<int>>& areas)
{
	ifstream in(filename);

	int n, m;
	in >> n >> m;
	areas.resize(n, vector<int>(m));

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			in >> areas[i][j];
		}
	}

	in.close();
}

int CountBreakPoints(vector<Interval>& intervals)
{
	int k = 0;
	for (auto& interval : intervals)
		k += interval.n;
	k++;

	return k;
}

void IntervalNumbering(vector<Interval>& intervals, int k)
{
	intervals[0].beginN = 0;
	int size = intervals.size();
	for (int i = 0; i < size - 1; i++)
	{
		intervals[i].endN = intervals[i].beginN + intervals[i].n * k;
		intervals[i + 1].beginN = intervals[i].endN;
	}

	intervals[size - 1].endN = intervals[size - 1].beginN + intervals[size - 1].n * k;
}

void BuildMesh(vector<Interval>& intervals, int k, vector<double>& x, vector<double>& hx)
{
	int pos = 0;
	for (auto interval : intervals)
	{
		double begin = interval.begin;
		double end = interval.end;
		int n = interval.n;

		x[pos++] = begin;

		double h = (end - begin) / n;
		for (int i = 1; i < n; i++, pos++)
		{
			double xi = begin + i * h;
			x[pos] = xi;
			hx[pos - 1] = xi - x[pos - 1];
		}
		hx[pos - 1] = end - x[pos - 1];
	}

	x[pos] = intervals.back().end;
	hx[pos - 1] = x[pos] - x[pos - 1];
}