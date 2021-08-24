using System;

namespace MathUtility
{
	public class VectorUtils
	{
		public static double Error(double[] a, double[] b)
		{
			double result = 0.0;
			int N = a.Length;

			for (int i = 0; i < N; i++)
				result += (a[i] - b[i]) * (a[i] - b[i]);

			return Math.Sqrt(result);
		}

		public static double DotProduct(double[] a, double[] b)
		{
			if (a.Length != b.Length)
				throw new Exception("vectors have different length");

			double result = 0.0;
			for (int i = 0; i < a.Length; i++)
				result += a[i] * b[i];

			return result;
		}

		public static double Norm(double[] a)
		{
			double result = 0.0;

			foreach (var value in a)
				result += value * value;

			return Math.Sqrt(result);
		}

		public static double RelativeError(double[] a, double[] b)
		{
			double result = 0.0;
			int N = a.Length;

			for (int i = 0; i < N; i++)
				result += (a[i] - b[i]) * (a[i] - b[i]);

			return Math.Sqrt(result) / Norm(b);
		}
	}
}
