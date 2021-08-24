using System;

namespace MathUtility
{
	public struct Point
	{
		public double X { get; set; }
		public double Y { get; set; }

		public Point(double x, double y)
		{
			X = x;
			Y = y;
		}

		public static Point operator +(Point a, Point b) => new Point(a.X + b.X, a.Y + b.Y);

		public static Point operator *(Point a, double constant) => new Point(a.X * constant, a.Y * constant);

		public static Point operator /(Point a, double constant) => new Point(a.X / constant, a.Y / constant);

		public static double Distance(Point a, Point b) => Math.Sqrt((a.X - b.X) * (a.X - b.X) + (a.Y - b.Y) * (a.Y - b.Y));

		public static Point Parse(string value)
		{
			string[] tokens = value.Split(' ');

			double.TryParse(tokens[0], out double a);
			double.TryParse(tokens[1], out double b);

			return new Point(a, b);
		}
	}
}
