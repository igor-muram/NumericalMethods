#pragma once

#include <fstream>
#include <string>
#include <vector>
#include <exception>
#include <cmath>

#include "MathUtility.h"

#pragma region SubArea

struct SubArea
{
	double x1, x2;
	double y1, y2;

	double mu;
	double J;
	int material;
};

std::istream& operator>>(std::istream& stream, SubArea& area)
{
	stream >> area.x1 >> area.x2 >> area.y1 >> area.y2 >> area.mu >> area.J >> area.material;
	return stream;
}

#pragma endregion

#pragma region Interval

struct Interval
{
	double a = 0.0;
	double  b = 0.0;
	double q = 0.0;
	int dir = 0;

	bool isUniform = false;

	size_t pointsCount = 0;
	double firstStep = 0.0;

	void Split()
	{
		if (MathUtility::IsEqual(a, b) || MathUtility::IsEqual(q, 0.0) || dir == 0)
			throw new std::exception("ERROR! Interval parameters are invalid");

		isUniform = q == 1.0;

		if (!isUniform)
		{
			int n = (int)(std::log((b - a) * (q - 1) / firstStep + 1) / std::log(q));

			double x1 = a + firstStep * (std::pow(q, n - 1) - 1) / (q - 1);
			double x2 = a + firstStep * (std::pow(q, n) - 1) / (q - 1);
			double lastStep = x2 - x1;


			if (b - x2 < lastStep * 0.5)
				pointsCount = n + 1;
			else
				pointsCount = n + 2;
		}
		else
		{
			int n = (int)((b - a) / firstStep);
			pointsCount = n + 1;
		}
	}

	double operator[] (int i)
	{
		if (i > pointsCount)
			throw new std::exception("ERROR! Invalid index");

		if (isUniform)
		{
			if (dir == 1)
				return GetUniformPoint(i);
			else
				return GetUniformPoint(pointsCount - i - 1);
		}
		else
		{
			if (dir == 1)
				return GetNonUniformPoint(i);
			else
				return GetNonUniformPoint(pointsCount - i - 1);
		}
	}

	double GetUniformPoint(int i)
	{
		return a + firstStep * i;
	}

	double GetNonUniformPoint(int i)
	{
		if (i != pointsCount - 1)
			return a + firstStep * (std::pow(q, i) - 1) / (q - 1);
		else
			return b;
	}
};

#pragma endregion

#pragma region Rectangle

struct Rectangle
{
	double x1, x2;
	double y1, y2;

	double mu;
	double J;
	int material;
};

#pragma endregion

using Mesh = std::vector<Rectangle>;

class Area
{
public:
	std::vector<SubArea> SubAreas;
	std::vector<Interval> XIntervals;
	std::vector<Interval> YIntervals;
	int XFactor = 0;
	int YFactor = 0;

public:
	void BuildMesh(Mesh& mesh)
	{
		if (SubAreas.size() == 0 || XIntervals.size() == 0 || YIntervals.size() == 0 || XFactor == 0 || YFactor == 0)
			throw new std::exception("ERROR! Area data is not initialized");

		for (auto& xInt : XIntervals)
			for (auto& yInt : YIntervals)
			{
				for (int i = 0; i < xInt.pointsCount - 1; i++)
				{
					for (int j = 0; j < yInt.pointsCount - 1; j++)
					{
						Rectangle rect;

						rect.x1 = xInt[i];
						rect.x2 = xInt[i + 1];
						rect.y1 = yInt[j];
						rect.y2 = yInt[j + 1];

						auto subarea = FindSubArea(xInt, yInt);

						rect.mu = subarea->mu;
						rect.J = subarea->J;
						rect.material = subarea->material;

						mesh.push_back(rect);
					}
				}
			}
	}

	void Clear()
	{
		SubAreas.clear();
		SubAreas.shrink_to_fit();

		XIntervals.clear();
		XIntervals.shrink_to_fit();

		YIntervals.clear();
		YIntervals.shrink_to_fit();
	}

private:
	SubArea* FindSubArea(Interval& x, Interval& y)
	{
		// Check if any subarea is the rect with X and Y sizes
		for (auto& sub : SubAreas)
		{
			if (MathUtility::IsEqual(x.a, sub.x1) && MathUtility::IsEqual(x.b, sub.x2))
				if (MathUtility::IsEqual(y.a, sub.y1) && MathUtility::IsEqual(y.b, sub.y2))
					return &sub;
		}

		// Check if any subarea contains rect with X and Y sizes
		for (auto& sub : SubAreas)
		{
			if (MathUtility::IsEqualOrLess(sub.x1, x.a) && MathUtility::IsEqualOrLess(x.b, sub.x2))
				if (MathUtility::IsEqualOrLess(sub.y1, y.a) && MathUtility::IsEqualOrLess(y.b, sub.y2))
					return &sub;
		}

		return nullptr;
	}

public:
	static Area& FromFile(std::string filename)
	{
		std::ifstream file(filename);

		if (file.is_open())
		{
			// Subareas reading
			int rectCount = 0;
			file >> rectCount;

			Area* area = new Area();

			for (int i = 0; i < rectCount; i++)
			{
				SubArea subarea;
				file >> subarea;

				area->SubAreas.push_back(subarea);
			}

			// Interval's info reading
			ReadIntervals(file, area->XIntervals);
			ReadIntervals(file, area->YIntervals);

			// splitting factors reading
			file >> area->XFactor >> area->YFactor;

			file.close();

			return *area;
		}
		else
		{
			std::string message = "ERROR! Failed to open file" + filename;
			throw new std::exception(message.c_str());
		}
	}

private:
	static void ReadIntervals(std::ifstream& file, std::vector<Interval>& intervals)
	{
		double initPoint = 0.0;
		int intervalNum = 0;
		file >> initPoint >> intervalNum;

		intervals.resize(intervalNum);

		for (int i = 0; i < 4; i++)
		{
			for (int k = 0; k < intervalNum; k++)
			{
				switch (i)
				{
				case 0:
				{
					double x = 0.0;
					file >> x;
					intervals[k].a = initPoint;
					intervals[k].b = x;
					initPoint = x;
					break;
				}
				case 1:
				{
					double step = 0.0;
					file >> step;
					intervals[k].firstStep = step;
					break;
				}
				case 2:
				{
					double q = 0.0;
					file >> q;
					intervals[k].q = q;
					break;
				}
				case 3:
				{
					int dir = 0;
					file >> dir;
					intervals[k].dir = dir;
					break;
				}

				default:
					break;
				}
			}
		}

	}
};