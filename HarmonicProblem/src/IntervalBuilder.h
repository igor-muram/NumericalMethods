#pragma once

#include "stdafx.h"
#include "FEMInfo.h"

class IntervalBuilder
{
public:
	IntervalBuilder(std::vector<RawInterval>& intervals)
	{
		for (auto i : intervals)
		{

			double begin = i.begin;
			double end = i.end;
			int n = i.n;
			double q = i.q;

			double h = (q == 1.0) ? ((end - begin) / n) : ((end - begin) * (1.0 - q) / (1.0 - pow(q, n)));

			double point = begin;
			for (int i = 0; i < n; i++)
			{
				points.push_back({ point, count++ });
				h *= q;
				point += h;
			}
		}

		points.push_back({ intervals.back().end, count });
		count++;
	}

	PI Begin() { return points.begin(); }
	PI End() { return points.end(); }

	int Count() { return count; }

	double operator[](int i)
	{
		return points[i].value;
	}

	double Back() { return points.back().value; }

private:
	std::vector<Point> points;
	int count = 0;
};