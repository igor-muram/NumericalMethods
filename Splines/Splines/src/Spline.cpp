#include "Spline.h"

void NaturalCubicSpline(vector<Point>& points, vector<Point>& spline)
{
	spline.clear();
	size_t n = points.size() - 1;

	vector<double> a(n + 1), b(n), d(n);

	for (size_t i = 0; i < n + 1; i++)
		a[i] = points[i].y;

	vector<double> h(n);
	for (size_t i = 0; i < n; i++)
		h[i] = points[i + 1].x - points[i].x;

	vector<double> alpha(n);
	for (size_t i = 1; i < n; i++)
		alpha[i] = (a[i + 1] - a[i]) * 3 / h[i] - (a[i] - a[i - 1]) * 3 / h[i - 1];

	vector<double> c(n + 1), l(n + 1), m(n + 1), z(n + 1);
	l[0] = 1.0;
	m[0] = z[0] = 0.0;

	for (size_t i = 1; i < n; i++)
	{
		l[i] = 2 * (points[i + 1].x - points[i - 1].x) - h[i - 1] * m[i - 1];
		m[i] = h[i] / l[i];
		z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
	}

	l[n] = 1;
	z[n] = c[n] = 0;

	for (int j = n - 1; j >= 0; j--)
	{
		c[j] = z[j] - m[j] * c[j + 1];
		b[j] = (a[j + 1] - a[j]) / h[j] - h[j] * (c[j + 1] + 2 * c[j]) / 3;
		d[j] = (c[j + 1] - c[j]) / (3 * h[j]);
	}

	for (size_t i = 0; i < n; i++)
	{
		double x0 = points[i].x;
		double x1 = points[i + 1].x;

		double ai = a[i];
		double bi = b[i];
		double ci = c[i];
		double di = d[i];
		double xi = points[i].x;

		function<double(double)> s = [=](double x) { return ai + bi * (x - xi) + ci * (x - xi) * (x - xi) + di * (x - xi) * (x - xi) * (x - xi); };
		double h = (x1 - x0) / 30;
		for (double x = x0; x < x1; x += h)
			spline.push_back(Point(x, s(x)));
	}

	double end = points.back().x;

	spline.push_back(
		Point(
			end,
			a[n - 1] +
			b[n - 1] * (end - points[n - 1].x) +
			c[n - 1] * (end - points[n - 1].x) * (end - points[n - 1].x) +
			d[n - 1] * (end - points[n - 1].x) * (end - points[n - 1].x) * (end - points[n - 1].x)));
}

void LagrangePolynomial(vector<Point>& points, vector<Point>& spline)
{
	spline.clear();
	size_t size = points.size();
	vector<function<double(double)>> basis;

	double x0 = 0;

	for (auto point : points)
	{
		x0 = point.x;
		function<double(double)> L = [&points, x0](double x)
		{
			double sum = 1.0;
			for (Point point : points)
				if (abs(x0 - point.x) > 1.0e-5)
					sum *= (x - point.x) / (x0 - point.x);

			return sum;
		};

		basis.push_back(L);
	}

	double a = points[0].x;
	double b = points.back().x;

	double h = (b - a) / 100;

	for (double x = a; x < b; x += h)
	{
		double s = 0;

		for (int i = 0; i < size; i++)
			s += basis[i](x) * points[i].y;

		spline.push_back(Point(x, s));
	}

	double s = 0;

	for (int i = 0; i < size; i++)
		s += basis[i](b) * points[i].y;

	spline.push_back(Point(b, s));
}

void LagrangeQuadratic(vector<Point>& points, vector<Point>& spline)
{
	spline.clear();

	size_t n = points.size() - 1;
	vector<double> h(n), beta(n);

	for (int i = 0; i < n; i++)
		h[i] = points[i + 1].x - points[i].x;

	vector<double> a(n), b(n), c(n);

	c[0] = (points[1].y - points[0].y) / (h[0] * h[0]) - beta[0] / h[0];
	b[0] = beta[0] - 2 * c[0] * points[0].x;
	a[0] = points[0].y - beta[0] * points[0].x + c[0] * points[0].x * points[0].x;

	for (int i = 1; i < n; i++)
	{
		beta[i] = b[i - 1] + 2 * c[i - 1] * points[i].x;
		c[i] = (points[i + 1].y - points[i].y) / (h[i] * h[i]) - beta[i] / h[i];
		b[i] = beta[i] - 2 * c[i] * points[i].x;
		a[i] = points[i].y - beta[i] * points[i].x + c[i] * points[i].x * points[i].x;
	}

	for (size_t i = 0; i < n; i++)
	{
		double x0 = points[i].x;
		double x1 = points[i + 1].x;

		double ai = a[i];
		double bi = b[i];
		double ci = c[i];

		function<double(double)> s = [=](double x) { return ai + bi * x + ci * x * x; };
		double h = (x1 - x0) / 30;
		for (double x = x0; x < x1; x += h)
			spline.push_back(Point(x, s(x)));
	}

	double end = points.back().x;
	spline.push_back(Point(end, a[n - 1] + b[n - 1] * end + c[n - 1] * end * end));
}

uint64_t Factorial(int n, int m = 1)
{
	if (n < 0)
		return 0;

	if (n == 0)
		return 1;

	uint64_t res = 1;
	for (int i = m; i <= n; i++)
		res *= i;

	return res;
}

uint64_t Binomial(int k, int n)
{
	if (k > n)
		return 0;

	uint64_t result = 0;
	int sub = n - k;

	if (k > sub)
	{
		result = Factorial(n, k + 1);
		result /= Factorial(sub);
	}
	else
	{
		result = Factorial(n, sub + 1);
		result /= Factorial(k);
	}

	return result;
}

void BezierCurve(vector<Point>& points, vector<Point>& spline)
{
	spline.clear();
	size_t n = points.size() - 1;
	uint32_t size = (n + 1) * 30;

	vector<double> coeffs;
	for (int i = 0; i <= n; i++)
		coeffs.push_back(Binomial(i, n));

	function<Point(double)> curve = [&points, &coeffs, n](double t)
	{
		Point C(0.0, 0.0);

		for (int k = 0; k <= n; k++)
			C = C + points[k] * coeffs[k] * pow((1 - t), n - k) * pow(t, k);

		return C;
	};

	double h = 1.0 / size;
	for (double t = 0.0; t < 1.0; t += h)
		spline.push_back(curve(t));
}

void Coordinates(int i, vector<Point>& points, vector<double>& coeffs)
{
	double x1 = points[i - 1].x;
	double x2 = points[i].x;
	double x3 = points[i + 1].x;
	double x4 = points[i + 2].x;

	double y1 = points[i - 1].y;
	double y2 = points[i].y;
	double y3 = points[i + 1].y;
	double y4 = points[i + 2].y;

	coeffs[0] = (-x1 + 3 * (x2 - x3) + x4) / 6.0;
	coeffs[1] = (-y1 + 3 * (y2 - y3) + y4) / 6.0;

	coeffs[2] = (x1 - 2 * x2 + x3) / 2.0;
	coeffs[3] = (y1 - 2 * y2 + y3) / 2.0;

	coeffs[4] = (x3 - x1) / 2.0;
	coeffs[5] = (y3 - y1) / 2.0;

	coeffs[6] = (x1 + 4 * x2 + x3) / 6.0;
	coeffs[7] = (y1 + 4 * y2 + y3) / 6.0;
}

void BSpline(vector<Point>& points, vector<Point>& spline)
{
	spline.clear();
	size_t n = points.size() - 1;

	vector<double> coeffs;
	coeffs.resize(8);

	int N = 100;

	if (n == 0)
		points.push_back(points[0]);

	if (n > 2)
	{
		spline.push_back(points[0]);
		points.push_back(points[n]);

		for (int i = 1; i < n; i++)
		{
			Coordinates(i, points, coeffs);
			for (int j = 0; j <= N; j++)
			{
				double t = (double)j / N;
				double X = ((coeffs[0] * t + coeffs[2]) * t + coeffs[4]) * t + coeffs[6];
				double Y = ((coeffs[1] * t + coeffs[3]) * t + coeffs[5]) * t + coeffs[7];

				spline.push_back(Point(X, Y));
			}
		}

		spline.push_back(points[n]);
		points.pop_back();
	}
}