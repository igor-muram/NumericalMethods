#pragma once

#include <vector>
#include <functional>
#include <cstdint>
#include <string>
#include <sstream>

using namespace std;

/// <summary>
/// Класс, представляющий точку
/// </summary>
struct Point
{
	Point() : x(0.0), y(0.0) {};
	Point(double x, double y) : x(x), y(y) {};
	double x, y;

	Point operator+(Point p)
	{
		return Point(p.x + x, p.y + y);
	}

	Point operator*(double a)
	{
		return Point(a * x, a * y);
	}

	string ToString()
	{
		stringstream stream;
		stream << "(" << x << ", " << y << ")";
		return stream.str();
	}
};

/// <summary>
/// Кусочно-кубическая интерполяция с помощью метода прогонки (натуральный кубический сплайн)
/// </summary>
/// <param name="points">Массив входных точек (x, y)</param>
/// <param name="spline">Массив точек для вывода сплайна</param>
void NaturalCubicSpline(vector<Point>& points, vector<Point>& spline);

/// <summary>
/// Интерполяционный многочлен Лагранжа
/// </summary>
/// <param name="points">Массив входных точек (x, y)</param>
/// <param name="spline">Массив точек для вывода сплайна</param>
void LagrangePolynomial(vector<Point>& points, vector<Point>& spline);

/// <summary>
/// Кусочно-квадратичная интерполяция
/// </summary>
/// <param name="points">Массив входных точек (x, y)</param>
/// <param name="spline">Массив точек для вывода сплайна</param>
void LagrangeQuadratic(vector<Point>& points, vector<Point>& spline);

/// <summary>
/// Кривая Безье
/// </summary>
/// <param name="points">Массив входных точек (x, y)</param>
/// <param name="spline">Массив точек для вывода сплайна</param>
void BezierCurve(vector<Point>& points, vector<Point>& spline);

/// <summary>
/// B-сплайн
/// </summary>
/// <param name="points">Массив входных точек (x, y)</param>
/// <param name="spline">Массив точек для вывода сплайна</param>
void BSpline(vector<Point>& points, vector<Point>& spline);