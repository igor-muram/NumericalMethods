#include "MeshBuilder.h"

void ReadIntervals(string filename, vector<Interval>& intervals)
{
	ifstream in(filename);

	int count;
	in >> count;

	for (int i = 0; i < count; i++)
	{
		Interval interval;
		in >> interval.begin >> interval.end >> interval.n >> interval.q;
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
		for (int j = 0; j < m; j++)
			in >> areas[i][j];

	in.close();
}

void ReadBoundaryConds(string filename, vector<BoundaryCondition>& conds)
{
	ifstream in(filename);

	int count;
	in >> count;
	conds.resize(count);

	for (int i = 0; i < count; i++)
	{
		in >> conds[i].xBegin >> conds[i].xEnd;
		in >> conds[i].yBegin >> conds[i].yEnd;
		in >> conds[i].functionNo;
	}

	in.close();
}

int CountNodes(vector<Interval>& intervals)
{
	int k = 0;
	for (auto& interval : intervals)
		k += interval.n;
	k++;

	return k;
}

void IntervalNumbering(vector<Interval>& intervals)
{
	intervals[0].beginN = 0;
	int size = intervals.size();
	for (int i = 0; i < size - 1; i++)
	{
		intervals[i].endN = intervals[i].beginN + intervals[i].n;
		intervals[i + 1].beginN = intervals[i].endN;
	}

	intervals[size - 1].endN = intervals[size - 1].beginN + intervals[size - 1].n;
}

void BoundaryCondsNumbering(vector<Interval>& intervalsX, vector<Interval>& intervalsY, vector<BoundaryCondition>& conds)
{
	vector<int> xMap, yMap;

	xMap.push_back(intervalsX[0].beginN);
	for (auto interval : intervalsX)
		xMap.push_back(interval.endN);

	yMap.push_back(intervalsY[0].beginN);
	for (auto interval : intervalsY)
		yMap.push_back(interval.endN);

	for (auto& cond : conds)
	{
		cond.xBegin = xMap[cond.xBegin];
		cond.yBegin = yMap[cond.yBegin];
		cond.xEnd = xMap[cond.xEnd];
		cond.yEnd = yMap[cond.yEnd];
	}
}

void BuildMesh(vector<Interval>& intervals, vector<double>& x)
{
	int pos = 0;
	for (auto interval : intervals)
	{
		double begin = interval.begin;
		double end = interval.end;
		int n = interval.n;
		double q = interval.q;

		double h = (q == 1.0) ? ((end - begin) / n) : ((end - begin) * (1.0 - q) / (1.0 - pow(q, n)));

		x[pos++] = begin;
		for (int i = 1; i < n; i++, pos++)
		{
			x[pos] = begin + i * h;
			h *= q;
		}
	}

	x[pos] = intervals.back().end;
}