using MathUtilities;
using SlaeSolver;
using System;
using System.Collections.Generic;
using FEM;
using System.Numerics;

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

		public static double[] DirectProblem(Source source, List<Receiver> receivers, double sigma1, double sigma2, double h1, double eps)
		{
			Dictionary<int, Material> materials = new Dictionary<int, Material>()
			{
				{ 0, new Material(sigma1) },
				{ 1, new Material(sigma2) }
			};

			Dictionary<AreaSide, (ConditionType, Func<double, double, double>)> Conditions = new Dictionary<AreaSide, (ConditionType, Func<double, double, double>)>()
			{
				{ AreaSide.TopFirst, (ConditionType.Second, (double r, double z) => 0) },
				{ AreaSide.Bottom,   (ConditionType.First, (double r, double z) => 0) },
				{ AreaSide.Right,    (ConditionType.First, (double r, double z) => 0) }
			};

			AreaInfoWithoutDelta areainfo = new AreaInfoWithoutDelta();
			areainfo.R0 = 0.0;
			areainfo.Z0 = 0.0;
			areainfo.Width = 10000;
			areainfo.FirstLayerHeight = h1;
			areainfo.SecondLayerHeight = 10000 - h1;

			areainfo.HorizontalStartStep = 0.1;
			areainfo.HorizontalCoefficient = 1.1;
			areainfo.VerticalStartStep = 0.1;
			areainfo.VerticalCoefficient = 1.1;
			areainfo.Materials = materials;
			areainfo.Conditions = Conditions;

			areainfo.SplitPoint = 0.5;

			GridBuilderWithoutDelta gb = new GridBuilderWithoutDelta(areainfo);
			gb.Build();

			foreach (Edge edge in gb.SB.Edges)
				edge.Function = (double r, double z) => 1.0 / (Math.PI * gb.ClosestSplitPoint * gb.ClosestSplitPoint);

			FEMProblemInfo info = new FEMProblemInfo();
			info.Points = gb.Points.ToArray();
			info.Materials = materials;
			info.Mesh = gb.Grid;
			info.FB = gb.FB;
			info.SB = gb.SB;
			info.F = (double r, double z) => 0.0;
			info.SolverType = SolverTypes.LOSLU;

			FEMrz fem = new FEMrz(info);
			fem.Solver.Eps = eps;
			fem.Solve();

			int n = receivers.Count;
			double[] V = new double[n];

			for (int i = 0; i < n; i++)
			{
				Vector3 M = receivers[i].M;
				Vector3 N = receivers[i].N;

				double vABM = fem.U(new Point(M.X - 100.0, M.Y)) - fem.U(new Point(M.X, M.Y));
				double vABN = fem.U(new Point(N.X - 100.0, N.Y)) - fem.U(new Point(N.X, N.Y));

				V[i] = vABM - vABN;
			}

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
