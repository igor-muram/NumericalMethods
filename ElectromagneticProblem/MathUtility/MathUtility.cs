using System;

namespace MathUtility
{
	public class FloatComparision
	{
		public static bool IsEqual(double a, double b)
		{
			return Math.Abs(b - a) < 1.0e-9;
		}

		public static bool IsEqualOrLess(double a, double b)
		{
			return Math.Abs(b - a) < 1.0e-9 || a < b;
		}

		public static bool IsEqualOrGreater(double a, double b)
		{
			return Math.Abs(b - a) < 1.0e-9 || a > b;
		}
	}
}
