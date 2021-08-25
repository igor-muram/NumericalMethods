#pragma once
#include "stdafx.h"
#include "FEMInfo.h"
#include "Matrix.h"

class PortraitBuilder
{
public:
	PortraitBuilder(int n, FEIterator begin, FEIterator end) : n(n), begin(begin), end(end)
	{
		connections.resize(n);
		BuildConnections();
	}

	void BuildSparse(SparseMatrix& A)
	{
		A.IA.resize(n + 1);
		A.IA[0] = A.IA[1] = 0;

		for (int i = 2; i < n + 1; i++)
			A.IA[i] = A.IA[i - 1] + connections[i - 1].size();

		for (int i = 0; i < n; i++)
		{
			for (auto j : connections[i])
			{
				A.JA.push_back(j);
			}
		}

		A.AL.resize(A.IA.back());
		A.AU.resize(A.IA.back());
	}

	void BuildProfile(ProfileMatrix& A)
	{
		A.IA.resize(n + 1);
		A.IA[0] = A.IA[1] = 0;

		for (int i = 2; i < n + 1; i++)
		{
			int width = i - *(connections[i - 1].begin()) - 1;
			A.IA[i] = A.IA[i - 1] + width;
		}


		A.AL.resize(A.IA.back());
		A.AU.resize(A.IA.back());
	}

private:
	int n, JASize;
	FEIterator begin, end;
	std::vector <std::set<int>> connections;

	void BuildConnections()
	{
		for (FEIterator iter = begin; iter != end; iter++)
		{
			FE e = *iter;
			for (auto v1: e.verts)
				for (auto v2 : e.verts)
				{
					if (v1 == v2)
					{
						connections[2 * v1 + 1].insert(2 * v1);
					}
					else
					{
						int a = v1;
						int b = v2;
						if (a > b) std::swap(a, b);

						connections[2 * b].insert(2 * a);
						connections[2 * b].insert(2 * a + 1);
						connections[2 * b + 1].insert(2 * a);
						connections[2 * b + 1].insert(2 * a + 1);
					}
				}
		}
	}
};
