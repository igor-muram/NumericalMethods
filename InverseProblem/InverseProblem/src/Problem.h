#pragma once

#include "Constants.h"

#include "Gauss.h"
#include "Classes.h"

#include "Matrix.h"
#include "Solvers.h"

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

double Difference(std::vector<Receiver>& g, std::vector<Receiver>& realG)
{
	int N = g.size();

	double result = 0.0;
	for (int i = 0; i < N; i++)
		result += (g[i].value - realG[i].value) * (g[i].value - realG[i].value);

	return result / N;
}