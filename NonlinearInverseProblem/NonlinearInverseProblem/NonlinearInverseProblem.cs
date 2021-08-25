using MathUtilities;
using SlaeSolver;
using System;
using System.Collections.Generic;

namespace NonlinearInverseProblem
{
	public class Problem
	{
		public struct ProblemInfo
		{
			public Source Source;
			public List<Receiver> Receivers;
			public double[] V;
			public double[] TrueV;
			public double[] p;
		}

		public static double[] DirectProblem(Source source, List<Receiver> receivers, double p)
		{
			double[] V = new double[receivers.Count];
			int N = receivers.Count;

			for (int i = 0; i < N; i++)
				V[i] = source.Potential(receivers[i], p);

			return V;
		}

		public static double[] InverseProblem(ProblemInfo info)
		{
			int N = info.p.Length;

			double[] dp = new double[N];
			FullSparseMatrix A = new FullSparseMatrix(N);

			for (int i = 0; i < N; i++)
			{
				for (int j = 0; j < N; j++)
				{
					double value = 0.0;
					int k = 0;
					foreach (var receiver in info.Receivers)
					{
						double w = 1.0 / info.TrueV[k];
						Func<double, double> func = (double p) => info.Source.Potential(receiver, p);
						double derivative = Utilities.Derivative(func, info.p[i * N + j], 0.0001);

						value += (w * derivative) * (w * derivative);
						k++;
					}

					A.Set(i, j, value);
				}
			}

			double[] B = new double[N];

			for (int i = 0; i < N; i++)
			{
				double value = 0.0;
				int k = 0;
				foreach (var receiver in info.Receivers)
				{
					double w = 1.0 / info.TrueV[k];

					Func<double, double> func = (double p) => info.Source.Potential(receiver, p);
					double derivative = Utilities.Derivative(func, info.p[i], 0.0001);
					value += w * w * derivative * (info.V[k] - info.TrueV[k]);
					k++;
				}

				B[i] = -value;
			}

			LOSLU solver = new LOSLU();
			return solver.Solve(A, B);
		}

		public static double Functional(double[] V, double[] TrueV)
		{
			double result = 0.0;
			int N = V.Length;

			for (int i = 0; i < N; i++)
			{
				double w = 1.0 / TrueV[i];
				result += w * (V[i] - TrueV[i]) * w * (V[i] - TrueV[i]);
			}

			return result;
		}
	}
}
