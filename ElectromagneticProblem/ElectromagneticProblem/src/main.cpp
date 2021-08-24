#include "Area.h"

#include <iostream>

int main()
{
	//Area area = Area::FromFile("input/grid.txt");

	//Mesh mesh = area.BuildMesh();

	Interval i;
	i.a = 0;
	i.b = 16;
	i.initStep = 1;
	i.q = 2;
	i.qDirection = 1;

	i.splitter = IntervalSplitter(i.a, i.b, i.initStep, i.q, i.qDirection);

	for (int k = 0; k < i.splitter.pointsCount; k++)
		std::cout << i.splitter[k] << "\t";

	return 0;
}