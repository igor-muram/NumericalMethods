#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <sstream>

#include "Problem.h"
#include "CSV.h"

const int ReceiversCount = 200;

const int ElementCount = 200;
const int ElementCountX = 20;
const int ElementCountZ = 10;

const double ElementWidth = 100.0;
const double ElementHeight = 50.0;

const double ReceiversX0 = -1000.0;
const double ReceiversZ0 = 0.0;
const double ReceiversStep = 10.0;

void CreateGrid(std::vector<Rect>& grid)
{
	for (int i = 0; i < ElementCountZ; i++)
		for (int j = 0; j < ElementCountX; j++)
			grid.emplace_back(ReceiversX0 + j * ElementWidth, ReceiversX0 + (j + 1) * ElementWidth, i * -ElementHeight, (i + 1) * -ElementHeight);
}

void CreateReceivers(std::vector<Receiver>& receivers)
{
	for (int i = 0; i < ReceiversCount; i++)
		receivers.emplace_back(ReceiversX0 + i * ReceiversStep, ReceiversZ0);
}

void CreateRealParameters(std::vector<Element>& real_Element, std::vector<Rect>& grid)
{
	for (int i = 0; i < ElementCount; i++)
		real_Element.emplace_back(grid[i], 0.0);

	real_Element[49].value = 2.0;
	real_Element[50].value = 1.0;
	real_Element[69].value = 1.0;
	real_Element[70].value = 2.0;

	CSV::PrintElements("RealValues.csv", real_Element, ElementCountX, ElementCountZ);
}

void CreateInitialParameters(std::vector<Element>& init_parameters, std::vector<Rect>& grid)
{
	for (int i = 0; i < ElementCount; i++)
		init_parameters.emplace_back(grid[i], 0.0);

	init_parameters[49].value = 1.0;
	init_parameters[50].value = 1.5;
	init_parameters[69].value = 1.5;
	init_parameters[70].value = 1.0;
}

void CreateGammaInfo(std::vector<std::pair<double, int>>& gammaInfo)
{
	gammaInfo.resize(ElementCount);

	gammaInfo[0] = std::make_pair(10.0e-15, 3);
	gammaInfo[ElementCountX - 1] = std::make_pair(10.0e-15, 3);
	gammaInfo[ElementCount - ElementCountX] = std::make_pair(10.0e-15, 3);
	gammaInfo[ElementCount - 1] = std::make_pair(10.0e-15, 3);

	for (int i = 1; i < ElementCountX - 1; i++)
	{
		gammaInfo[i] = std::make_pair(10.0e-15, 5);
		gammaInfo[(ElementCountZ - 1) * ElementCountX + i] = std::make_pair(10.0e-15, 5);
	}

	for (int j = ElementCountX ; j < (ElementCountX - 2) * ElementCountZ; j += ElementCountX)
	{
		gammaInfo[j] = std::make_pair(10.0e-15, 5);
		gammaInfo[j + ElementCountX - 1] = std::make_pair(10.0e-15, 5);
	}

	for (int i = 1; i < ElementCountZ - 1; i++)
		for (int j = 1; j < ElementCountX - 1; j++)
			gammaInfo[i * ElementCountX + j] = std::make_pair(10.0e-15, 8);
}

void TestWithSmoothingRegularization()
{
	// Grid
	std::vector<Rect> grid;
	CreateGrid(grid);

	// True p
	std::vector<Element> real_Element;
	CreateRealParameters(real_Element, grid);

	// True g
	std::vector<Receiver> real_g;
	CreateReceivers(real_g);
	DirectProblem(real_Element, real_g);

	// Initial p
	std::vector<Element> init_Element;
	CreateInitialParameters(init_Element, grid);

	// Initital g
	std::vector<Receiver> init_g;
	CreateReceivers(init_g);
	DirectProblem(init_Element, init_g);

	// Initial gamma
	std::vector<std::pair<double, int>> gammaInfo;
	CreateGammaInfo(gammaInfo);

	std::vector<Element> changed_elements;
	CreateInitialParameters(changed_elements, grid);

	double alpha = 1.0e-17;
	bool badConds = false;
	int iter = 0;

	ProblemInfo info;
	info.alpha = 1.0e-17;
	info.g0 = &init_g;
	info.realG = &real_g;
	info.grid = &grid;

	std::vector<double> dp;
	RegularizedInverseProblem(init_g, real_g, grid, dp, 1.0e-17);
	// Get new parameters
	for (int i = 0; i < ElementCount; i++)
		changed_elements[i].value = init_Element[i].value + dp[i];

	CSV::PrintElements("Result (Without Smoothing).csv", changed_elements, ElementCountX, ElementCountZ);

	do
	{
		// Create G matrix
		Matrix G;
		CreateGMatrix(gammaInfo, ElementCountX, G);
		info.G = &G;

		// Calculate delta p
		SmoothingRegularizedInverseProblem(info, dp);

		// Get new parameters
		for (int i = 0; i < ElementCount; i++)
			changed_elements[i].value = init_Element[i].value + dp[i];

		std::stringstream stream;
		stream << "Result (" << iter++ << ").csv";

		CSV::PrintElements(stream.str(), changed_elements, ElementCountX, ElementCountZ);

		// bad smoothing conditions
		std::set<std::pair<int, int>> badCondsCoords;

		// Check smoothing
		badConds = CheckSmoothingCondition(changed_elements, badCondsCoords, ElementCountX, 0.1);

		ChangeGammas(gammaInfo, badCondsCoords);
	} while (badConds);
}

int main()
{
	TestWithSmoothingRegularization();
	return 0;
}