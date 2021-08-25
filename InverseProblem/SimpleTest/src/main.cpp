#include <iostream>
#include <vector>

#include "Problem.h"

const int ReceiversCount = 2;
const int ElementsCount = 4;
const double ReceiversX0 = -5.0;
const double ReceiversY0 = 0.0;
const double ReceiversStep = 10.0;

void CreateGrid(std::vector<Rect>& grid)
{
	grid.reserve(4);
	grid.emplace_back(-100, 0, -150, -100);
	grid.emplace_back(-100, 0, -200, -150);
	grid.emplace_back(0, 100, -150, -100);
	grid.emplace_back(0, 100, -200, -150);
}

void CreateRealParameters(std::vector<Element>& real_elements, std::vector<Rect>& grid)
{
	real_elements.reserve(ElementsCount);
	real_elements.emplace_back(grid[0], 2);
	real_elements.emplace_back(grid[1], 1);
	real_elements.emplace_back(grid[2], 1);
	real_elements.emplace_back(grid[3], 2);
}

void CreateInitialParameters(std::vector<Element>& init_parameters, std::vector<Rect>& grid)
{
	init_parameters.reserve(ElementsCount);
	init_parameters.emplace_back(grid[0], 1.1);
	init_parameters.emplace_back(grid[1], 1.5);
	init_parameters.emplace_back(grid[2], 1.5);
	init_parameters.emplace_back(grid[3], 1.1);
}

void TestWithoutRegularization()
{
	// Grid
	std::vector<Rect> grid;
	grid.reserve(4);
	grid.emplace_back(-100, 0, -150, -100);
	grid.emplace_back(-100, 0, -200, -150);
	grid.emplace_back(0, 100, -150, -100);
	grid.emplace_back(0, 100, -200, -150);


	// True p
	std::vector<Element> real_elements;
	real_elements.reserve(4);
	real_elements.emplace_back(grid[0], 2);
	real_elements.emplace_back(grid[1], 1);
	real_elements.emplace_back(grid[2], 1);
	real_elements.emplace_back(grid[3], 2);

	// True g
	std::vector<Receiver> real_g;
	for (int i = 0; i < ReceiversCount; i++)
		real_g.emplace_back(-1000.0 + i * 10.0, 0.0);

	DirectProblem(real_elements, real_g);

	// Initial p
	std::vector<Element> elements0;
	elements0.reserve(4);
	elements0.emplace_back(grid[0], 10.0);
	elements0.emplace_back(grid[1], 10.0);
	elements0.emplace_back(grid[2], 10.0);
	elements0.emplace_back(grid[3], 10.0);

	// Initital g
	std::vector<Receiver> g0;
	for (int i = 0; i < ReceiversCount; i++)
		g0.emplace_back(-1000.0 + i * 10.0, 0.0);
	DirectProblem(elements0, g0);


	// Delta p
	std::vector<double> dp;
	InverseProblem(g0, real_g, grid, dp);

	// Change p0
	for (int i = 0; i < ElementsCount; i++)
		elements0[i].value += dp[i];

	// g1
	std::vector<Receiver> g1;
	for (int i = 0; i < ReceiversCount; i++)
		g1.emplace_back(-1000.0 + i * 10.0, 0.0);
	DirectProblem(elements0, g1);

	std::cout << "Difference: " << Difference(g1, real_g) << std::endl;
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
	for (int i = 0; i < ReceiversCount; i++)
		real_g.emplace_back(ReceiversX0 + i * ReceiversStep, ReceiversY0);

	DirectProblem(real_elements, real_g);

	// Initial p
	std::vector<Element> init_elements;
	CreateInitialParameters(init_elements, grid);

	// Initital g
	std::vector<Receiver> init_g;
	for (int i = 0; i < ReceiversCount; i++)
		init_g.emplace_back(ReceiversX0 + i * ReceiversStep, ReceiversY0);
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
		for (int i = 0; i < ReceiversCount; i++)
			new_g.emplace_back(ReceiversX0 + i * ReceiversStep, ReceiversY0);
		DirectProblem(changed_elements, new_g);

		// Calculate difference
		std::cout << "Difference: " << Difference(new_g, real_g) << "\t";
		std::cout << "Alpha: " << alpha << "\t";
		for (int i = 0; i < ElementsCount; i++)
			std::cout << changed_elements[i].value << "\t\t";

		std::cout << std::endl;

		// Change alpha
		alpha *= 10.0;
	}
}

int main()
{
	TestWithRegularization();
	return 0;
}