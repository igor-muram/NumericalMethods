#include <stdio.h>
#include <stdlib.h>
#include <math.h>

float F(float t, float y)
{
	return -25 * y + cos(t) + 25 * sin(t);
}

float f(float t)
{
	return sin(t)+ exp(-25 * t);
}

void ExplicitEuler(float a, float b, float h, float y0)
{
	int n = (b - a) / h;

	printf_s("%f\t\t%f\t\t%f\t\t%f\n", 0, y0, f(0), fabs(y0 - f(0)));

	for (int i = 0; i < n; i++)
	{
		float t = i * h;
		y0 += h * F(t, y0);

		printf_s("%f\t\t%f\t\t%f\t\t%f\n", t + h, y0, f(t + h), fabs(y0 - f(t + h)));
	}
}

void ImplicitEuler(float a, float b, float h, float y0, float eps)
{
	int n = (b - a) / h;
	printf_s("%f\t\t%f\t\t%f\t\t%f\n", 0.0, y0, f(0), fabs(y0 - f(0)));

	for (int i = 0; i < n; i++)
	{
		float t = (i + 1) * h;
		float yk0 = y0, yk1 = yk0;
		do
		{
			yk0 = yk1;
			yk1 = yk0 - (y0 + h * F(t + h, yk0) - yk0) / (-25 * h - 1);
		} 
		while (fabs(yk1 - yk0) >= eps);
		y0 = yk1;

		printf_s("%f\t\t%f\t\t%f\t\t%f\n", t, y0, f(t), fabs(y0 - f(t)));
	}
}

void ModifiedEuler(float a, float b, float h, float y0)
{
	int n = (b - a) / h;

	printf_s("%f\t\t%f\t\t%f\t\t%f\n", 0.0, y0, f(0.0), fabs(y0 - f(0.0)));

	for (int i = 0; i < n; i++)
	{
		float t = i * h;
		y0 += h * (F(t, y0) + F(t + h, y0 + h * F(t, y0))) / 2;

		printf_s("%f\t\t%f\t\t%f\t\t%f\n", t + h, y0, f(t + h), fabs(y0 - f(t + h)));
	}
}

void Trapeze(float a, float b, float h, float y0, float eps)
{
	int n = (b - a) / h;
	printf_s("%f\t\t%f\t\t%f\t\t%f\n", 0.0, y0, f(0), fabs(y0 - f(0)));

	for (int i = 0; i < n; i++)
	{
		float t = (i + 1) * h;
		
		float yk0 = y0, yk1 = yk0;
		do
		{
			yk0 = yk1;
			yk1 = yk0 - (y0 + h * (F(t, y0) + F(t + h, yk0)) / 2 - yk0) / (-25 * h / 2 - 1);
		} 
		while (fabs(yk1 - yk0) >= eps);

		y0 = yk1;

		printf_s("%f\t\t%f\t\t%f\t\t%f\n", t, y0, f(t), fabs(y0 - f(t)));
	}
}

void Runge(float a, float b, float h, float y0)
{
	int n = (b - a) / h;

	for (int i = 0; i < n; i++)
	{
		float t = i * h;
		float k1 = F(t, y0);
		float k2 = F(t + h / 2, y0 + h * k1 / 2);
		float k3 = F(t + h / 2, y0 + h * k2 / 2);
		float k4 = F(t + h, y0 + h * k3);

		y0 += h * (k1 + 2 * k2 + 2 * k3 + k4) / 6;

		printf_s("%f\t\t%f\t\t%f\t\t%f\n", t, y0, f(t + h), fabs(y0 - f(t + h)));
	}
}

int main()
{
	Trapeze(0.0, 2.0, 0.2, 1.0, 1.0e-4);
	system("pause");
	return 0;
}