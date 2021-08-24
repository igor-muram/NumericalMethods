#include "Area.h"

#include <iostream>

int main()
{
	//Area area = Area::FromFile("input/grid.txt");

	//Mesh mesh = area.BuildMesh();

	Interval i;
	i.a = 0;
	i.b = 9.5;
	i.firstStep = 1;
	i.q = 2;
	i.dir = 1;
	i.Split();

	for (int k = 0; k < i.pointsCount; k++)
		std::cout << i[k] << "\t";

	return 0;
}