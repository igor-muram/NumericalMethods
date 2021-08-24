using System;

namespace MathUtility
{
	public class Optimization
	{
		public static double GoldenRatio(Func<double, double> f, double a, double b, double eps = 1.0e-7)
		{
			int k = 0;
			double diff = (b - a);

			int n = (int)(Math.Log((b - a) / eps) / Math.Log((Math.Sqrt(5) + 1) / 2));
			double prevDiff, ratio;

			double x1 = a + 0.381966011 * diff;
			double x2 = b - 0.381966011 * diff;

			double f1 = f(x1);
			double f2 = f(x2);
			k += 2;

			for (int i = 0; i <= n; i++)
			{
				if (f1 < f2)
				{
					b = x2;
					x2 = x1;
					x1 = a + (b - a) * 0.381966011;
					f2 = f1;
					f1 = f(x1);
					k++;
				}
				else
				{
					a = x1;
					x1 = x2;
					x2 = b - (b - a) * 0.381966011;
					f1 = f2;
					f2 = f(x2);
					k++;
				}

				prevDiff = diff;
				diff = (b - a);

				ratio = prevDiff / diff;
			}

			return (x1 + x2) / 2;
		}
	}
}
