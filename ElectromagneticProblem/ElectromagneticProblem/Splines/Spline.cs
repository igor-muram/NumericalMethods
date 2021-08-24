using System;
using System.Globalization;
using System.IO;
using System.Linq;

namespace ElectromagneticProblem.Splines
{
	public interface ISpline
	{
		public double GetValue(double arg);
		public double GetDerivative(double arg);
	}

	public class CubicHermiteSpline : ISpline
	{
		int count;
		double[] args;
		double[] values;
		double[] derivatives;

		CubicHermiteSpline() { }

		void CalculateDerivatives()
		{
			double h1 = args[1] - args[0];
			double h2 = args[2] - args[1];

			double l = -(2 * h1 + h2) / (h1 * (h1 + h2));
			double c = (h1 + h2) / (h1 * h2);
			double r = -h1 / (h2 * (h1 + h2));

			derivatives[0] = values[0] * l + values[1] * c + values[2] * r;

			h1 = args[count - 2] - args[count - 3];
			h2 = args[count - 1] - args[count - 2];

			l = h1 / (h1 * (h1 + h2));
			c = -(h1 + h2) / (h1 * h2);
			r = (2 * h2 + h1) / (h2 * (h1 + h2));

			derivatives[count - 1] = values[count - 3] * l + values[count - 2] * c + values[count - 1] * r;

			for (int i = 1; i < count - 1; i++)
			{
				h1 = args[i] - args[i - 1];
				h2 = args[i + 1] - args[i];

				l = -h2 / (h1 * (h1 + h2));
				c = -(h1 - h2) / (h1 * h2);
				r = h1 / (h2 * (h1 + h2));
				derivatives[i] = values[i - 1] * l + values[i] * c + values[i + 1] * r;
			}
		}

		static public CubicHermiteSpline FromFile(string filename)
		{
			string[] data;
			try
			{
				data = File.ReadAllLines(filename);
			}
			catch
			{
				Console.WriteLine($"ERROR! Failed to open file {filename}");
				return null;
			}

			CubicHermiteSpline spline = new CubicHermiteSpline();

			int count = int.Parse(data[0], NumberStyles.Any);

			spline.count = count;
			spline.args = new double[count];
			spline.values = new double[count];
			spline.derivatives = new double[count];

			for (int i = 0; i < count; i++)
			{
				string[] tokens = data[i + 1].Split(' ', StringSplitOptions.RemoveEmptyEntries);

				double arg = double.Parse(tokens[1], NumberStyles.Any, CultureInfo.InvariantCulture);
				double value = double.Parse(tokens[0], NumberStyles.Any, CultureInfo.InvariantCulture);

				spline.args[i] = arg;
				spline.values[i] = value;
			}

			spline.CalculateDerivatives();

			return spline;
		}

		public double GetValue(double arg)
		{
			if (arg > args[count - 1] || arg < args[0])
				throw new Exception("ERROR! Bad arg for spline");

			int i1 = 0;
			int i2 = 0;

			for (int i = 1; i < count; i++)
				if (arg >= args[i - 1] && arg <= args[i])
				{
					i1 = i - 1;
					i2 = i;
					break;
				}

			double h = args[i2] - args[i1];
			double ksi = (arg - args[i1]) / h;

			double phi0 = 1 - 3 * ksi * ksi + 2 * ksi * ksi * ksi;
			double phi1 = h * (ksi - 2 * ksi * ksi + ksi * ksi * ksi);
			double phi2 = 3 * ksi * ksi - 2 * ksi * ksi * ksi;
			double phi3 = h * (-ksi * ksi + ksi * ksi * ksi);

			double value = values[i1] * phi0 + values[i2] * phi2 + derivatives[i1] * phi1 + derivatives[i2] * phi3;
			return value;
		}

		public double GetDerivative(double arg)
		{
			if (arg > args[count - 1] || arg < args[0])
				throw new Exception("ERROR! Bad arg for spline");

			int i1 = 0;
			int i2 = 0;

			for (int i = 1; i < count; i++)
				if (arg >= args[i - 1] && arg <= args[i])
				{
					i1 = i - 1;
					i2 = i;
					break;
				}

			double h = args[i2] - args[i1];
			double ksi = (arg - args[i1]) / h;

			double phi0 = h * (-6 * ksi + 6 * ksi * ksi);
			double phi1 = h * h * (1 - 4 * ksi + 3 * ksi * ksi);
			double phi2 = h * (6 * ksi - 6 * ksi * ksi);
			double phi3 = h * h * (-2 * ksi + 3 * ksi * ksi);

			double value = values[i1] * phi0 + values[i2] * phi2 + derivatives[i1] * phi1 + derivatives[i2] * phi3;
			return value;
		}
	}

	public class MuSpline : ISpline
	{
		int count;
		double[] args;
		double[] values;
		double[] derivatives;

		MuSpline() { }

		void CalculateDerivatives()
		{
			double h1 = args[1] - args[0];
			double h2 = args[2] - args[1];

			double l = -(2 * h1 + h2) / (h1 * (h1 + h2));
			double c = (h1 + h2) / (h1 * h2);
			double r = -h1 / (h2 * (h1 + h2));

			derivatives[0] = values[0] * l + values[1] * c + values[2] * r;

			h1 = args[count - 2] - args[count - 3];
			h2 = args[count - 1] - args[count - 2];

			l = h1 / (h1 * (h1 + h2));
			c = -(h1 + h2) / (h1 * h2);
			r = (2 * h2 + h1) / (h2 * (h1 + h2));

			derivatives[count - 1] = values[count - 3] * l + values[count - 2] * c + values[count - 1] * r;

			for (int i = 1; i < count - 1; i++)
			{
				h1 = args[i] - args[i - 1];
				h2 = args[i + 1] - args[i];

				l = -h2 / (h1 * (h1 + h2));
				c = -(h1 - h2) / (h1 * h2);
				r = h1 / (h2 * (h1 + h2));
				derivatives[i] = values[i - 1] * l + values[i] * c + values[i + 1] * r;
			}
		}

		static public MuSpline FromFile(string filename)
		{
			string[] data;
			try
			{
				data = File.ReadAllLines(filename);
			}
			catch
			{
				Console.WriteLine($"ERROR! Failed to open file {filename}");
				return null;
			}

			MuSpline spline = new MuSpline();

			int count = int.Parse(data[0], NumberStyles.Any);

			spline.count = count;
			spline.args = new double[count];
			spline.values = new double[count];
			spline.derivatives = new double[count];

			for (int i = 0; i < count; i++)
			{
				string[] tokens = data[i + 1].Split(' ', StringSplitOptions.RemoveEmptyEntries);

				double arg = double.Parse(tokens[1], NumberStyles.Any, CultureInfo.InvariantCulture);
				double value = double.Parse(tokens[0], NumberStyles.Any, CultureInfo.InvariantCulture);

				spline.args[i] = arg;
				spline.values[i] = value;
			}

			spline.CalculateDerivatives();

			return spline;
		}

		public double GetValue(double arg)
		{
			if (arg < args[0])
				throw new Exception("ERROR! Bad arg for spline");

			if (arg > args[count - 1])
			{
				double B = arg;
				double BN = args.Last();
				double muN = values.Last();

				double value = BN / B * (1.0 / muN - 1.0) + 1.0;
				return 1.0 / value;
			}

			int i1 = 0;
			int i2 = 0;

			for (int i = 1; i < count; i++)
				if (arg >= args[i - 1] && arg <= args[i])
				{
					i1 = i - 1;
					i2 = i;
					break;
				}

			double h = args[i2] - args[i1];
			double ksi = (arg - args[i1]) / h;

			double phi0 = 1 - 3 * ksi * ksi + 2 * ksi * ksi * ksi;
			double phi1 = h * (ksi - 2 * ksi * ksi + ksi * ksi * ksi);
			double phi2 = 3 * ksi * ksi - 2 * ksi * ksi * ksi;
			double phi3 = h * (-ksi * ksi + ksi * ksi * ksi);
	
			return values[i1] * phi0 + values[i2] * phi2 + derivatives[i1] * phi1 + derivatives[i2] * phi3;
		}

		public double GetDerivative(double arg)
		{
			if (arg < args[0])
				throw new Exception("ERROR! Bad arg for spline");

			if (arg > args[count - 1])
			{
				double B = arg;
				double BN = args.Last();
				double muN = values.Last();

				double value1 = BN / B * (1.0 / muN - 1.0) + 1.0;
				value1 *= value1;

				double value2 = BN / (B * B) * (1.0 / muN - 1.0);

				return value2 / value1;
			}

			int i1 = 0;
			int i2 = 0;

			for (int i = 1; i < count; i++)
				if (arg >= args[i - 1] && arg <= args[i])
				{
					i1 = i - 1;
					i2 = i;
					break;
				}

			double h = args[i2] - args[i1];
			double ksi = (arg - args[i1]) / h;

			double phi0 = h * (-6 * ksi + 6 * ksi * ksi);
			double phi1 = h * h * (1 - 4 * ksi + 3 * ksi * ksi);
			double phi2 = h * (6 * ksi - 6 * ksi * ksi);
			double phi3 = h * h * (-2 * ksi + 3 * ksi * ksi);
		
			return values[i1] * phi0 + values[i2] * phi2 + derivatives[i1] * phi1 + derivatives[i2] * phi3;
		}
	}
}
