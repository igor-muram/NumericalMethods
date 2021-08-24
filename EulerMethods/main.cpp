#include <stdio.h>
#include <stdlib.h>
#include <math.h>

float F(float t, float y)
{
	return 2 * t * y;
}

float f(float t)
{
	return exp(t * t);
}

void ImplicitEuler(float a, float b, float h, float y0, float eps)
{
	int n = (b - a) / h;
	printf_s("%f\t\t%f\t\t%f\t\t%f\n", 0.0, y0, f(0), fabs(y0 - f(0)));

	for (int i = 0; i < n; i++)
	{
		float t = (i + 1) * h;
		float yk0 = y0;
		float yk1 = y0 + h * F(t, yk0);
		while (fabs(yk0 - yk1) >= eps)
		{
			yk0 = yk1;
			yk1 = y0 + h * F(t, yk0);
		}

		y0 = yk1;
		printf_s("%f\t\t%f\t\t%f\t\t%f\n", t, y0, f(t), fabs(y0 - f(t)));
	}
}

void Trapeze(float a, float b, float h, float y0, float eps)
{
	int n = (b - a) / h;
	printf_s("%f\t\t%f\t\t%f\t\t%f\n", 0.0, y0, f(0), fabs(y0 - f(0)));

	for (int i = 0; i < n; i++)
	{
		float t = (i + 1) * h;
		float yk0 = y0;
		float yk1 = y0 + h * (F(t, yk0) + F(t - h, y0)) / 2;
		while (fabs(yk0 - yk1) >= eps)
		{
			yk0 = yk1;
			yk1 = y0 + h * (F(t, yk0) + F(t - h, y0)) / 2;
		}

		y0 = yk1;
		printf_s("%f\t\t%f\t\t%f\t\t%f\n", t, y0, f(t), fabs(y0 - f(t)));
	}
}

int main()
{
	ImplicitEuler(0, 1, 0.000001, 1.0, 10e-12);
	printf_s("\n\n\n");

	system("pause");
	return 0;
}