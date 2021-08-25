#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <sstream>

#include "Problem.h"
#include "CSV.h"

const int ReceiversCount = 200;

const int ElementsCount = 200;
const int ElementsCountX = 20;
const int ElementsCountZ = 10;

const double ElementWidth = 100.0;
const double ElementHeight = 50.0;

const double ReceiversX0 = -1000.0;
const double ReceiversZ0 = 0.0;
const double ReceiversStep = 10.0;

void CreateGrid(std::vector<Rect>& grid)
{
	for (int i = 0; i < ElementsCountZ; i++)
		for (int j = 0; j < ElementsCountX; j++)
			grid.emplace_back(ReceiversX0 + j * ElementWidth, ReceiversX0 + (j + 1) * ElementWidth, i * -ElementHeight, (i + 1) * -ElementHeight);
}

void CreateReceivers(std::vector<Receiver>& receivers)
{
	for (int i = 0; i < ReceiversCount; i++)
		receivers.emplace_back(ReceiversX0 + i * ReceiversStep, ReceiversZ0);
}

void CreateRealParameters(std::vector<Element>& real_elements, std::vector<Rect>& grid)
{
	for (int i = 0; i < ElementsCount; i++)
		real_elements.emplace_back(grid[i], 0.0);

	real_elements[49].value = 2.0;
	real_elements[50].value = 1.0;
	real_elements[69].value = 1.0;
	real_elements[70].value = 2.0;

	CSV::PrintElements("RealValues.csv", real_elements, ElementsCountX, ElementsCountZ);
}

void CreateInitialParameters(std::vector<Element>& init_parameters, std::vector<Rect>& grid)
{
	for (int i = 0; i < ElementsCount; i++)
		init_parameters.emplace_back(grid[i], 0.0);

	init_parameters[49].value = 1.0;
	init_parameters[50].value = 1.5;
	init_parameters[69].value = 1.5;
	init_parameters[70].value = 1.0;
}

void TestWithRegularization()
{
	// Grid
	std::vector<Rect> grid;
	CreateGrid(grid);

	// True p
	std::vector<Element> real_elements;
	CreateRealParameters(real_elements, grid);

	// True g
	std::vector<Receiver> real_g;
	CreateReceivers(real_g);
	DirectProblem(real_elements, real_g);

	// Initial p
	std::vector<Element> init_elements;
	CreateInitialParameters(init_elements, grid);

	// Initital g
	std::vector<Receiver> init_g;
	CreateReceivers(init_g);
	DirectProblem(init_elements, init_g);

	std::vector<Element> changed_elements;
	CreateInitialParameters(changed_elements, grid);

	
	double alpha = 1.0e-20;
	for (int i = 0; i < 15; i++)
	{
		// Calculate delta p
		std::vector<double> dp;
		RegularizedInverseProblem(init_g, real_g, grid, dp, alpha);

		// Get new parameters
		for (int i = 0; i < ElementsCount; i++)
			changed_elements[i].value = init_elements[i].value + dp[i];

		// Calculate new g
		std::vector<Receiver> new_g;
		CreateReceivers(new_g);
		DirectProblem(changed_elements, new_g);

		// Calculate difference
		std::cout << "Difference: " << Difference(new_g, real_g) << "\t";
		std::cout << "Alpha: " << alpha << "\t";

		std::cout << setprecision(4) <<
			changed_elements[49].value << "\t" <<
			changed_elements[50].value << "\t" <<
			changed_elements[69].value << "\t" <<
			changed_elements[70].value << std::endl;

		std::cout << std::endl;

		std::stringstream stream;
		stream << "Result (" << alpha << ").csv";
		CSV::PrintElements(stream.str(), { changed_elements, ElementsCountX, ElementsCountZ, alpha,  Difference(new_g, real_g) });


		// Change alpha
		alpha *= 10.0;
	}
}

int main()
{
	TestWithRegularization();
	return 0;
}