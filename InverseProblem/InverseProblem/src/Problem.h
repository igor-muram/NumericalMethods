#pragma once

#include "Constants.h"

#include "Gauss.h"
#include "Classes.h"

#include "Matrix.h"
#include "Solvers.h"

#include <set>
#include <map>

void DirectProblem(std::vector<Element>& elements, std::vector<Receiver>& receivers)
{
	size_t size = receivers.size();
	for (Receiver& receiver : receivers)
	{
		receiver.value = 0.0;
		for (auto& element : elements)
		{
			std::function<double(double, double)> f = [&receiver](double x, double z)
			{
				double recx = receiver.x;
				double recz = receiver.z;
				double r = sqrt((x - recx) * (x - recx) + (z - recz) * (z - recz));
				return z / (r * r * r);
			};

			receiver.value += element.value * Gauss2(f, element.rect) / (4 * PI);
		}
	}
}

void InverseProblem(std::vector<Receiver>& g0, std::vector<Receiver>& realG, std::vector<Rect>& Grid, std::vector<double>& dp)
{
	size_t gridSize = Grid.size();
	size_t receiversCount = g0.size();
	Matrix A(gridSize);

	for (int q = 0; q < gridSize; q++)
		for (int s = 0; s < gridSize; s++)
		{
			for (Receiver& receiver : g0)
			{
				std::function<double(double, double)> f = [&receiver](double x, double z)
				{
					double recx = receiver.x;
					double recz = receiver.z;
					double r = sqrt((x - recx) * (x - recx) + (z - recz) * (z - recz));
					return z / (r * r * r);
				};

				A(q, s) += Gauss2(f, Grid[q]) * Gauss2(f, Grid[s]) / (16.0 * PI * PI);
			}
		}

	vector<double> b(gridSize, 0.0);

	for (int q = 0; q < gridSize; q++)
	{
		for (int i = 0; i < receiversCount; i++)
		{
			std::function<double(double, double)> f = [i, &g0](double x, double z)
			{
				double recx = g0[i].x;
				double recz = g0[i].z;
				double r = sqrt((x - recx) * (x - recx) + (z - recz) * (z - recz));
				return z / (r * r * r);
			};

			b[q] += (g0[i].value - realG[i].value) * Gauss2(f, Grid[q]) / (4 * PI);
		}

		b[q] *= -1;
	}

	dp = std::vector<double>(gridSize, 0.0);

	Solvers::LOS(A, dp, b);

}

void RegularizedInverseProblem(std::vector<Receiver>& g0, std::vector<Receiver>& realG, std::vector<Rect>& Grid, std::vector<double>& dp, double alpha)
{
	size_t gridSize = Grid.size();
	size_t receiversCount = g0.size();
	Matrix A(gridSize);

	for (int q = 0; q < gridSize; q++)
		for (int s = 0; s < gridSize; s++)
		{
			for (Receiver& receiver : g0)
			{
				std::function<double(double, double)> f = [&receiver](double x, double z)
				{
					double recx = receiver.x;
					double recz = receiver.z;
					double r = sqrt((x - recx) * (x - recx) + (z - recz) * (z - recz));
					return z / (r * r * r);
				};

				A(q, s) += Gauss2(f, Grid[q]) * Gauss2(f, Grid[s]) / (16.0 * PI * PI);
			}
		}

	for (int i = 0; i < gridSize; i++)
		A.DI[i] += alpha;

	vector<double> b(gridSize, 0.0);

	for (int q = 0; q < gridSize; q++)
	{
		for (int i = 0; i < receiversCount; i++)
		{
			std::function<double(double, double)> f = [i, &g0](double x, double z)
			{
				double recx = g0[i].x;
				double recz = g0[i].z;
				double r = sqrt((x - recx) * (x - recx) + (z - recz) * (z - recz));
				return z / (r * r * r);
			};

			b[q] += (g0[i].value - realG[i].value) * Gauss2(f, Grid[q]) / (4 * PI);
		}

		b[q] *= -1;
	}

	dp = std::vector<double>(gridSize, 0.0);

	Solvers::LOS(A, dp, b);

}

struct ProblemInfo
{
	std::vector<Receiver>* g0 = nullptr;
	std::vector<Receiver>* realG = nullptr;
	std::vector<Rect>* grid = nullptr;
	double alpha = 0.0;
	Matrix* G;
};

void SmoothingRegularizedInverseProblem(ProblemInfo& info, std::vector<double>& dp)
{
	size_t gridSize = info.grid->size();
	size_t receiversCount = info.g0->size();
	Matrix A(gridSize);

	for (int q = 0; q < gridSize; q++)
		for (int s = 0; s < gridSize; s++)
		{
			for (Receiver& receiver : *info.g0)
			{
				std::function<double(double, double)> f = [&receiver](double x, double z)
				{
					double recx = receiver.x;
					double recz = receiver.z;
					double r = sqrt((x - recx) * (x - recx) + (z - recz) * (z - recz));
					return z / (r * r * r);
				};

				A(q, s) += Gauss2(f, info.grid->at(q)) * Gauss2(f, info.grid->at(s)) / (16.0 * PI * PI);
			}
		}

	for (int i = 0; i < gridSize; i++)
		A.DI[i] += info.alpha;

	for (int i = 0; i < gridSize; i++)
		for (int j = 0; j < gridSize; j++)
			A(i, j) += (*info.G)(i, j);

	vector<double> b(gridSize, 0.0);

	for (int q = 0; q < gridSize; q++)
	{
		for (int i = 0; i < receiversCount; i++)
		{
			std::function<double(double, double)> f = [i, info](double x, double z)
			{
				double recx = info.g0->at(i).x;
				double recz = info.g0->at(i).z;
				double r = sqrt((x - recx) * (x - recx) + (z - recz) * (z - recz));
				return z / (r * r * r);
			};

			b[q] += (info.g0->at(i).value - info.realG->at(i).value) * Gauss2(f, info.grid->at(q)) / (4 * PI);
		}

		b[q] *= -1;
	}

	dp = std::vector<double>(gridSize, 0.0);

	Solvers::LOS(A, dp, b);
}

void CreateGMatrix(std::vector<std::pair<double, int>>& gammaInfo, int ElementCountX, Matrix& G)
{
	size_t size = gammaInfo.size();
	G = Matrix(size);

	for (int i = 0; i < size; i++)
	{
		double sum = 0.0;
		for (int j = -1; j < 2; j++)
			for (int k = -1; k < 2; k++)
			{
				if (j != 0 || k != 0)
				{
					int index = i + j * ElementCountX + k;

					if (index > 0 && index < size)
						sum += gammaInfo[index].first;
				}
			}

		G(i, i) = gammaInfo[i].second * gammaInfo[i].first + sum;
	}

	for (int i = 0; i < size; i++)
		for (int j = 0; j < size; j++)
			if (i != j)
				G(i, j) = -(gammaInfo[i].first + gammaInfo[j].first);
}

bool CheckSmoothingCondition(std::vector<Element>& p, std::set<std::pair<int, int>>& badCondsCoords, int ElementCountX, double threshold)
{
	bool badCondsExists = false;

	size_t size = p.size();

	for (int i = 0; i < size; i++)
	{
		double value1 = p[i].value;

		if (abs(value1) > threshold)
		{
			for (int j = -1; j < 2; j++)
				for (int k = -1; k < 2; k++)
				{
					if (j != 0 || k != 0)
					{
						int near = i + j * ElementCountX + k;

						if (near >= 0 && near < size)
						{
							double value2 = p[near].value;

							if (value1 * value2 > 0 && abs(value2) > threshold)
							{
								double ratio = value1 / value2;
								if (ratio > 2.0 || ratio < 0.5)
									if (badCondsCoords.find(std::make_pair(near, i)) == badCondsCoords.end())
										badCondsCoords.insert(std::make_pair(i, near));
							}
						}
					}
				}
		}

	}

	if (badCondsCoords.size() > 0)
		badCondsExists = true;

	return badCondsExists;
}

void ChangeGammas(std::vector<std::pair<double, int>>& gammaInfo, std::set<std::pair<int, int>>& badCondsCoords)
{
	//std::map<int, bool> alreadyChanged;

	for (auto& pair : badCondsCoords)
	{
		int i = pair.first;
		int j = pair.second;
		
		gammaInfo[i].first *= 10.0;
		gammaInfo[j].first *= 10.0;
	}
}

double Difference(std::vector<Receiver>& g, std::vector<Receiver>& realG)
{
	int N = g.size();

	double result = 0.0;
	for (int i = 0; i < N; i++)
		result += (g[i].value - realG[i].value) * (g[i].value - realG[i].value);

	return result / N;
}