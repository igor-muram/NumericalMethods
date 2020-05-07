#pragma once

#include <vector>
#include "FEMInfo.h"
#include "Interval.h"

class MeshBuilder
{
public:
	MeshBuilder(std::vector<Interval>& intervals)
	{
		for (auto interval : intervals)
		{
			double begin = interval.begin;
			double end = interval.end;
			int n = interval.n;
			double q = interval.q;

			double h = (q == 1.0) ? ((end - begin) / n) : ((end - begin) * (1.0 - q) / (1.0 - pow(q, n)));

			FiniteElement e;
			e.begin = begin;
			e.end = begin + h;
			e.v.push_back(nodeCount++);
			e.v.push_back(nodeCount++);
			e.v.push_back(nodeCount);
			elements.push_back(e);

			for (int i = 1; i < n; i++)
			{
				h *= q;

				FiniteElement e;
				e.begin = elements.back().end;
				e.end = elements.back().end + h;
				e.v.push_back(nodeCount++);
				e.v.push_back(nodeCount++);
				e.v.push_back(nodeCount);
				elements.push_back(e);
			}
		}
		nodeCount++;
	}

	FEIterator Begin() { return elements.begin(); }
	FEIterator End() { return elements.end(); }

	int GetNodeCount() { return nodeCount; }

private:
	int nodeCount = 0;
	std::vector<FiniteElement> elements;
};