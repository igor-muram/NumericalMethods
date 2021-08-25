using System;
using System.Collections.Generic;

namespace MathUtilities
{
	public class Utilities
	{
		public static double Derivative(Func<double, double> func, double point, double delta)
		{
			return (func(point + delta) - func(point)) / delta;
		}

		public static double DotProduct(double[] a, double[] b)
		{
			if (a.Length != b.Length)
				throw new Exception("vectors have different length");

			double result = 0.0;
			int N = a.Length;

			for (int i = 0; i < N; i++)
				result += a[i] * b[i];

			return result;
		}

		public static int Factorial(int n)
		{
			int result = 1;
			for (int i = 2; i <= n; i++)
				result *= i;
			return result;
		}

		public static double Det(Point a, Point b, Point c)
		{
			double r1 = a.R;
			double z1 = a.Z;
			double r2 = b.R;
			double z2 = b.Z;
			double r3 = c.R;
			double z3 = c.Z;

			return (r2 - r1) * (z3 - z1) - (r3 - r1) * (z2 - z1);
		}

		public static double[] Alpha(Point a, Point b, Point c)
		{
			double[] alpha = new double[6];
			double r1 = a.R;
			double z1 = a.Z;
			double r2 = b.R;
			double z2 = b.Z;
			double r3 = c.R;
			double z3 = c.Z;

			double D = Det(a, b, c);

			alpha[0] = (z2 - z3) / D;
			alpha[1] = (r3 - r2) / D;
			alpha[2] = (z3 - z1) / D;
			alpha[3] = (r1 - r3) / D;
			alpha[4] = (z1 - z2) / D;
			alpha[5] = (r2 - r1) / D;

			return alpha;
		}

		public static Point[] TriangleCubicPoints(Point a, Point b, Point c)
		{
			Point[] points = new Point[10];

			points[0] = a;
			points[1] = b;
			points[2] = c;
			points[3] = (((b / 2.0) + a) * 2.0 / 3.0);
			points[4] = (((b * 2.0) + a) / 3.0);
			points[5] = (((c / 2.0) + a) * 2.0 / 3.0);
			points[6] = (((c * 2.0) + a) / 3.0);
			points[7] = (((c / 2.0) + b) * 2.0 / 3.0);
			points[8] = (((c * 2.0) + b) / 3.0);
			points[9] = ((a + b + c) / 3.0);

			return points;
		}

		public static double NewtonCotes(double a, double b, Func<double, double> f)
		{
			double h = (b - a) / 1000;
			double result = f(a);

			for (double x = h; x < b; x += h)
				result += 2 * f(x);

			result += f(b);
			result *= 7;

			for (double x = a; x < b; x += h)
				result += 32 * f(x + h / 4) + 12 * f(x + h / 2) + 32 * f(x + 3 * h / 4);

			result = result * 0.5 * h / 45;

			return result;
		}

		public static double[] Gauss(double[][] A, double[] b)
		{
			int N = A.Length;
			double[] x = new double[N];

			for (int k = 0; k < N - 1; k++)
			{
				// Поиск ведущего элемента
				double max = Math.Abs(A[k][k]);
				int m = k;
				for (int i = k + 1; i < N; i++)
					if (Math.Abs(A[k][i]) > max)
					{
						max = Math.Abs(A[k][i]);
						m = i;
					}

				// Обмен местами b[m] и b[k]
				(b[m], b[k]) = (b[k], b[m]);
				// Обмен местами k-ого и m-ого столбцов
				for (int j = k; j < N; j++)
					(A[k][j], A[m][j]) = (A[m][j], A[k][j]);

				// Обнуление k-ого столбца
				for (int i = k + 1; i < N; i++)
				{
					double t = A[i][k] / A[k][k];
					b[i] -= t * b[k];
					for (int j = k + 1; j < N; j++)
						A[i][j] -= t * A[k][j];
				}
			}

			// Вычисление вектора x
			x[N - 1] = b[N - 1] / A[N - 1][N - 1];
			for (int k = N - 2; k >= 0; k--)
			{
				double sum = 0;
				for (int j = k + 1; j < N; j++)
					sum += A[k][j] * x[j];

				x[k] = (b[k] - sum) / A[k][k];
			}

			return x;
		}
	}
}
