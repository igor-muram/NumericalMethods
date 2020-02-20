#pragma once

#include <vector>
#include <string>
#include <fstream>

using namespace std;

struct Interval
{
	double begin, end;
	int beginN, endN;
	int n;
};

struct Area
{
	Interval intervalX, intervalY;
};

const double eps = 1.0e-5;

void ReadIntervals(string filename, vector<Interval>& intervals);
void ReadAreaMatrix(string filename, vector<vector<int>>& areas);
int CountBreakPoints(vector<Interval>& intervals);
void IntervalNumbering(vector<Interval>& intervals, int k);
void BuildMesh(vector<Interval>& intervals, int k, vector<double>& x, vector<double>& h);