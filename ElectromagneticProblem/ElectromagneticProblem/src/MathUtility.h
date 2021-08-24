#pragma once

#include <cmath>

class MathUtility
{
public:
	static bool IsEqual(double a, double b)
	{
		return std::abs(b - a) < 1.0e-9;
	}

	static bool IsEqualOrLess(double a, double b)
	{
		return std::abs(b - a) < 1.0e-9 || a < b;
	}
};