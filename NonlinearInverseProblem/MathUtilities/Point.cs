using System;

namespace MathUtilities
{
	public struct Point
	{
		public double R { get; set; }
		public double Z { get; set; }

		public Point(double r, double z)
		{
			R = r;
			Z = z;
		}

		public static Point operator+(Point a, Point b)
		{
			return new Point(a.R + b.R, a.Z + b.Z);
		}

		public static Point operator*(Point a, double constant)
		{
			return new Point(a.R * constant, a.Z * constant);
		}

		public static Point operator/(Point a, double constant)
		{
			return new Point(a.R / constant, a.Z / constant);
		}

		public static double Distance(Point a, Point b)
		{
			return Math.Sqrt((a.R - b.R) * (a.R - b.R) + (a.Z - b.Z) * (a.Z - b.Z));
		}

		public static Point Parse(string value)
		{
			string[] tokens = value.Split(' ');

			double.TryParse(tokens[0], out double a);
			double.TryParse(tokens[1], out double b);

			return new Point(a, b);
		}
	}
}
