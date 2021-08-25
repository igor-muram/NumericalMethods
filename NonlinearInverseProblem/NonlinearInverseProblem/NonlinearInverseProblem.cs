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
			public double h;
			public double sigma1;
			public double sigma2;
		}

		public static FEMrz FEM(double sigma1, double sigma2, double h1, double eps)
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

			areainfo.HorizontalStartStep = 0.2;
			areainfo.HorizontalCoefficient = 1.1;
			areainfo.VerticalStartStep = 0.2;
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

			return fem;
		}

		public static double[] FEMDerivative(ProblemInfo info, double dh, double eps)
		{
			Console.WriteLine($"Start to calculate derivatives for {info.h}");

			int n = info.Receivers.Count;

			FEMrz fem1 = FEM(info.sigma1, info.sigma2, info.h, eps);
			FEMrz fem2 = FEM(info.sigma1, info.sigma2, info.h + dh, eps);

			double[] V1 = new double[n];
			for (int i = 0; i < n; i++)
			{
				Vector3 M = info.Receivers[i].M;
				Vector3 N = info.Receivers[i].N;

				double vABM = fem1.U(new Point(M.X - 100.0, M.Y)) - fem1.U(new Point(M.X, M.Y));
				double vABN = fem1.U(new Point(N.X - 100.0, N.Y)) - fem1.U(new Point(N.X, N.Y));

				V1[i] = vABM - vABN;
			}

			double[] V2 = new double[n];
			for (int i = 0; i < n; i++)
			{
				Vector3 M = info.Receivers[i].M;
				Vector3 N = info.Receivers[i].N;

				double vABM = fem2.U(new Point(M.X - 100.0, M.Y)) - fem2.U(new Point(M.X, M.Y));
				double vABN = fem2.U(new Point(N.X - 100.0, N.Y)) - fem2.U(new Point(N.X, N.Y));

				V2[i] = vABM - vABN;
			}

			double[] dV = new double[3];
			for (int i = 0; i < n; i++)
				dV[i] = (V2[i] - V1[i]) / dh;


			Console.WriteLine("Finished to calculate derivatives");
			Console.WriteLine();
			return dV;
		}

		public static double[] DirectProblem(ProblemInfo info, double eps)
		{
			FEMrz fem = FEM(info.sigma1, info.sigma2, info.h, eps);

			int n = info.Receivers.Count;
			double[] V = new double[n];

			for (int i = 0; i < n; i++)
			{
				Vector3 M = info.Receivers[i].M;
				Vector3 N = info.Receivers[i].N;

				double vABM = fem.U(new Point(M.X - 100.0, M.Y)) - fem.U(new Point(M.X, M.Y));
				double vABN = fem.U(new Point(N.X - 100.0, N.Y)) - fem.U(new Point(N.X, N.Y));

				V[i] = vABM - vABN;
			}

			return V;
		}

		public static double InverseProblem(ProblemInfo info, double eps)
		{
			double[] derivatives = FEMDerivative(info, 0.0001, eps);

			double A = 0.0;

			int k = 0;
			foreach (var receiver in info.Receivers)
			{
				double w = 1.0 / info.TrueV[k];
				A += (w * derivatives[k]) * (w * derivatives[k]);
				k++;
			}

			double B = 0.0;
			k = 0;
			foreach (var receiver in info.Receivers)
			{
				double w = 1.0 / info.TrueV[k];
				B -= w * w * derivatives[k] * (info.V[k] - info.TrueV[k]);
				k++;
			}

			return B / A;
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
