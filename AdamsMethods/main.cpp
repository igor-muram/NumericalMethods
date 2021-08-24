#include <cstdio>
#include <Functional>
#include <vector>

using namespace std;

double F(double t, double y)
{
	return 2 * t * y;
}

double f(double t)
{
	return exp(t * t);
}

double Runge(double a, double h, double y0)
{
	double k1 = F(a, y0);
	double k2 = F(a + h / 2, y0 + h * k1 / 2);
	double k3 = F(a + h / 2, y0 + h * k2 / 2);
	double k4 = F(a + h, y0 + h * k3);

	return y0 + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6;
}

void ExplicitAdams3Order(double a, double b, double h, double y0)
{
	int n = (b - a) / h;

	// Первые несколько значений расчитываем с помощью Рунге-Кутта
	vector<double> y0s;
	y0s.push_back(y0);
	for (int i = 0; i < 2; i++)
	{
		y0 = Runge(i * h, h, y0);
		y0s.push_back(y0);
	}

	// Вывод первых нескольких значений
	for (int i = 0; i < y0s.size(); i++)
		printf_s("%f\t\t%.7f\t\t%.7f\t\t%.7f\n", i * h, y0s[i], f(i * h), fabs(y0s[i] - f(i * h)));

	// Расчет остальных значений
	double y = y0s[2];
	for (int i = 0; i < n - 2; i++)
	{
		double t = (i + 2) * h;
		double temp = 
			23 * F(t, y0s[2]) -
			16 * F(t - h, y0s[1]) +
			5 * F(t - 2 * h, y0s[0]);

		y += h * temp / 12;

		y0s.insert(y0s.end(), y);
		y0s.erase(y0s.begin());

		printf_s("%f\t\t%.7f\t\t%.7f\t\t%.7f\n", t + h, y0s[2], f(t + h), fabs(y0s[2] - f(t + h)));
	}
}

void ExplicitAdams4Order(double a, double b, double h, double y0)
{
	int n = (b - a) / h;

	// Первые несколько значений расчитываем с помощью Рунге-Кутта
	vector<double> y0s;
	y0s.push_back(y0);
	for (int i = 0; i < 3; i++)
	{
		y0 = Runge(i * h, h, y0);
		y0s.push_back(y0);
	}

	// Вывод первых нескольких значений
	for (int i = 0; i < y0s.size(); i++)
		printf_s("%f\t\t%.7f\t\t%.7f\t\t%.7f\n", i * h, y0s[i], f(i * h), fabs(y0s[i] - f(i * h)));

	// Расчет остальных значений
	double y = y0s[3];
	for (int i = 0; i < n - 3; i++)
	{
		double t = (i + 3) * h;
		double temp =
			55 * F(t, y0s[3]) -
			59 * F(t - h, y0s[2]) +
			37 * F(t - 2 * h, y0s[1]) -
			9 * F(t - 3 * h, y0s[0]);

		y += h * temp / 24;

		y0s.insert(y0s.end(), y);
		y0s.erase(y0s.begin());

		printf_s("%f\t\t%.7f\t\t%.7f\t\t%.7f\n", t + h, y0s[3], f(t + h), fabs(y0s[3] - f(t + h)));
	}
}

void ImplicitAdams3Order(double a, double b, double h, double y0, double eps)
{
	int n = (b - a) / h;

	// Первые несколько значений расчитываем с помощью Рунге-Кутта
	vector<double> y0s;
	y0s.push_back(y0);
	y0s.push_back(Runge(0.0, h, y0));

	// Вывод первых нескольких значений
	for (int i = 0; i < y0s.size(); i++)
		printf_s("%f\t\t%.7f\t\t%.7f\t\t%.7f\n", i * h, y0s[i], f(i * h), fabs(y0s[i] - f(i * h)));

	// Расчет остальных значений
	for (int i = 0; i < n - 1; i++)
	{
		double t = (i + 1) * h;

		double temp =
			8 * F(t, y0s[1]) -
			F(t - h, y0s[0]);

		// Метод простой итерации для получения решения нелинейного уравнения
		double yk1 = y0s[1], yk0 = 0;
		while (fabs(yk1 - yk0) >= eps)
		{
			yk0 = yk1;
			yk1 = y0s[1] + h * (5 * F(t + h, yk0) + temp) / 12;
		}

		y0s.insert(y0s.end(), yk1);
		y0s.erase(y0s.begin());

		printf_s("%f\t\t%.7f\t\t%.7f\t\t%.7f\n", t + h, y0s[1], f(t + h), fabs(y0s[1] - f(t + h)));
	}
}

void ImplicitAdams4Order(double a, double b, double h, double y0, double eps)
{
	int n = (b - a) / h;

	// Первые несколько значений расчитываем с помощью Рунге-Кутта
	vector<double> y0s;
	y0s.push_back(y0);
	for (int i = 0; i < 2; i++)
	{
		y0 = Runge(i * h, h, y0);
		y0s.push_back(y0);
	}

	// Вывод первых нескольких значений
	for (int i = 0; i < y0s.size(); i++)
		printf_s("%f\t\t%.7f\t\t%.7f\t\t%.7f\n", i * h, y0s[i], f(i * h), fabs(y0s[i] - f(i * h)));

	// Расчет остальных значений
	for (int i = 0; i < n - 2; i++)
	{
		double t = (i + 2) * h;

		double temp =
			19 * F(t, y0s[2]) -
			5 * F(t - h, y0s[1]) +
			F(t - 2 * h, y0s[0]);

		// Метод простой итерации для получения решения нелинейного уравнения
		double yk1 = y0s[2], yk0 = 0;
		while (fabs(yk1 - yk0) >= eps)
		{
			yk0 = yk1;
			yk1 = y0s[2] + h * (9 * F(t + h, yk0) + temp) / 24;
		}

		y0s.insert(y0s.end(), yk1);
		y0s.erase(y0s.begin());

		printf_s("%f\t\t%.7f\t\t%.7f\t\t%.7f\n", t + h, y0s[2], f(t + h), fabs(y0s[2] - f(t + h)));
	}
}

int main()
{
	double a = 0.0;
	double b = 1.0;
	double h = 0.05;
	double y0 = 1.0;
	double eps = 1.0e-7;

	ImplicitAdams3(a, b, h, y0, eps);
	return 0;
}