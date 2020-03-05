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
	double q;
};

struct BoundaryCondition
{
	int xBegin, xEnd;
	int yBegin, yEnd;
	int functionNo;
};

void ReadIntervals(string filename, vector<Interval>& intervals);		
void ReadAreaMatrix(string filename, vector<vector<int>>& areas);	
void ReadBoundaryConds(string filename, vector<BoundaryCondition>& conds);
int CountNodes(vector<Interval>& intervals);							
void IntervalNumbering(vector<Interval>& intervals);					
void BoundaryCondsNumbering(vector<Interval>& intervalsX, vector<Interval>& intervalsY, vector<BoundaryCondition>& conds);
void BuildMesh(vector<Interval>& intervals, vector<double>& x);