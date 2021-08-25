#pragma once
#include "stdafx.h"
#include "FEMInfo.h"

class MeshBuilder
{
public:
	MeshBuilder(
		PI beginX, PI endX, int xCount, 
		PI beginY, PI endY, int yCount, 
		PI beginZ, PI endZ, int zCount)
	{
		connections.resize(xCount * yCount * zCount, -1);

		for (PI z = beginZ; z < endZ - 1; z++)
			for (PI y = beginY; y < endY - 1; y++)
				for (PI x = beginX; x < endX - 1; x++)
				{
					FE e;
					e.x[0] = x->value;
					e.y[0] = y->value;
					e.z[0] = z->value;
					e.x[1] = (x + 1)->value;
					e.y[1] = (y + 1)->value;
					e.z[1] = (z + 1)->value;

					int xs[2] = { x->i, x->i + 1 };
					int ys[2] = { y->i, y->i + 1 };
					int zs[2] = { z->i, z->i + 1 };

					for (int i = 0; i < 2; i++)
						for (int j = 0; j < 2; j++)
							for (int k = 0; k < 2; k++)
								e.verts.push_back(zs[i] * yCount * xCount + ys[j] * xCount + xs[k]);

					elements.push_back(e);
				}
			
	}

	FEIterator Begin() { return elements.begin(); }
	FEIterator End() { return elements.end(); }

private:
	std::vector<FE> elements;
	std::vector<int> connections;
	int nodeCount = 0;
};