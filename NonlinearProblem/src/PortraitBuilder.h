#pragma once

#include <vector>

#include "FEMInfo.h"
#include "Matrix.h"

class PortraitBuilder
{
public:
	PortraitBuilder(int nodeCount, FEIterator begin, FEIterator end): nodeCount(nodeCount), begin(begin), end(end)
	{
		connections.resize(nodeCount);
		BuildConnections();
	}

	void Build(Matrix& A)
	{
		A.IA.resize(nodeCount + 1);
		A.IA[0] = A.IA[1] = 0;

		for (int i = 2; i < nodeCount + 1; i++)
			A.IA[i] = A.IA[i - 1] + connections[i - 1].size();
	}

private:
	int nodeCount, JASize;
	FEIterator begin, end;
	std::vector <std::vector<int>> connections;

	void BuildConnections()
	{
		for (FEIterator iter = begin; iter != end; iter++)
		{
			FiniteElement e = *iter;
			connections[e.v[2]].push_back(e.v[0]);
			connections[e.v[2]].push_back(e.v[1]);
			connections[e.v[1]].push_back(e.v[0]);
		}
	}
};