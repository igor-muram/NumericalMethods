using SlaeSolver;
using FEM;
using System.Collections.Generic;
using System;
using System.Numerics;
using NonlinearInverseProblem;
using MathUtilities;

namespace ExampleTest
{
	class Program
	{
		static double sigma = 1.0;

		static void FEMTest()
		{
			Dictionary<int, Material> materials = new Dictionary<int, Material>()
			{
				{ 0, new Material(sigma) },
				{ 1, new Material(sigma) }
			};

			Dictionary<AreaSide, (ConditionType, Func<double, double, double>)> conditions = new Dictionary<AreaSide, (ConditionType, Func<double, double, double>)>()
			{
				{ AreaSide.Top,		(ConditionType.First, (double r, double z) => r * r * z) },
				{ AreaSide.Left,		(ConditionType.First, (double r, double z) => r * r * z) },
				{ AreaSide.Bottom,	(ConditionType.First, (double r, double z) => r * r * z) },
				{ AreaSide.Right,		(ConditionType.First, (double r, double z) => r * r * z) }
			};

			ImprovedAreaInfo areainfo = new ImprovedAreaInfo();
			areainfo.R0 = 0.0;
			areainfo.Z0 = 0.0;
			areainfo.Width = 100;
			areainfo.FirstLayerHeight = 100;
			areainfo.SecondLayerHeight = 300;

			areainfo.HorizontalStartStep = 10;
			areainfo.HorizontalCoefficient = 1.2;
			areainfo.VerticalStartStep = 10;
			areainfo.VerticalCoefficient = 1.2;
			areainfo.Materials = materials;
			areainfo.Conditions = conditions;

			ImprovedGridBuilder gb = new ImprovedGridBuilder(areainfo);
			gb.Build();

			Func<double, double, double> F = (double r, double z) => -sigma * 4.0 * z;

			ProblemInfo info = new ProblemInfo();
			info.Points = gb.Points.ToArray();
			info.Materials = materials;
			info.Mesh = gb.Grid;
			info.FB = gb.FB;
			info.F = F;
			info.SolverType = SolverTypes.LOSLU;
			FEMrz fem = new FEMrz(info);
			fem.Solve();

			Console.WriteLine(fem.Info.Points[35].R);
			Console.WriteLine(fem.Info.Points[35].Z);
			Console.WriteLine(fem.Weights[35]);
		}

		static void FEMDeltaTest()
		{
			Dictionary<int, Material> materials = new Dictionary<int, Material>()
			{
				{ 0, new Material(sigma) },
				{ 1, new Material(sigma) }
			};

			Dictionary<AreaSide, (ConditionType, Func<double, double, double>)> Conditions = new Dictionary<AreaSide, (ConditionType, Func<double, double, double>)>()
			{
				{ AreaSide.Top,		(ConditionType.SecondNull, (double r, double z) => 0) },
				{ AreaSide.Left,		(ConditionType.SecondNull, (double r, double z) => 0) },
				{ AreaSide.Bottom,	(ConditionType.First, (double r, double z) => 0) },
				{ AreaSide.Right,		(ConditionType.First, (double r, double z) => 0) }
			};

			ImprovedAreaInfo areainfo = new ImprovedAreaInfo();
			areainfo.R0 = 0.0;
			areainfo.Z0 = 0.0;
			areainfo.Width = 10000;
			areainfo.FirstLayerHeight = 100;
			areainfo.SecondLayerHeight = 9900;

			areainfo.HorizontalStartStep = 12.0;
			areainfo.HorizontalCoefficient = 1.0;
			areainfo.VerticalStartStep = 12.0;
			areainfo.VerticalCoefficient = 1.0;
			areainfo.Materials = materials;
			areainfo.Conditions = Conditions;

			ImprovedGridBuilder gb = new ImprovedGridBuilder(areainfo);
			gb.Build();

			ProblemInfoDelta infoDelta = new ProblemInfoDelta();
			infoDelta.Points = gb.Points.ToArray();
			infoDelta.Materials = materials;
			infoDelta.Mesh = gb.Grid;
			infoDelta.FB = gb.FB;
			infoDelta.R0 = 0.0;
			infoDelta.Z0 = 0.0;
			infoDelta.DeltaCoefficient = 1.0 / (2 * Math.PI);
			infoDelta.SolverType = SolverTypes.LOSLU;

			FEMrzDelta femDelta = new FEMrzDelta(infoDelta);
			femDelta.Solver.Eps = 1.0e-11;
			femDelta.Solve();

			Console.WriteLine();
			Console.WriteLine(femDelta.U(new Point(1000, 0)));
			Console.WriteLine();

			Console.WriteLine(femDelta.Solver.Difference);
			Console.WriteLine(femDelta.Solver.IterCount);
		}

		static void InverseProblemTest()
		{
			Vector3 A = new Vector3(0.0f, 0.0f, 0.0f);
			Vector3 B = new Vector3(100.0f, 0.0f, 0.0f);
			double I = 1.0;

			Source source = new Source(A, B, I);

			Vector3 M1 = new Vector3(200.0f, 0.0f, 0.0f);
			Vector3 N1 = new Vector3(300.0f, 0.0f, 0.0f);

			Vector3 M2 = new Vector3(500.0f, 0.0f, 0.0f);
			Vector3 N2 = new Vector3(600.0f, 0.0f, 0.0f);

			Vector3 M3 = new Vector3(1000.0f, 0.0f, 0.0f);
			Vector3 N3 = new Vector3(1100.0f, 0.0f, 0.0f);

			List<Receiver> receivers = new List<Receiver>() { new Receiver(M1, N1), new Receiver(M2, N2), new Receiver(M3, N3) };

			Problem.ProblemInfo info = new Problem.ProblemInfo();

			double trueP = 0.1;
			info.TrueV = Problem.DirectProblem(source, receivers, trueP);

			double initialP = 0.01;
			info.V = Problem.DirectProblem(source, receivers, initialP);

			info.p = new double[1] { initialP };
			info.Source = source;
			info.Receivers = receivers;

			while (Problem.Functional(info.V, info.TrueV) > 1.0e-7)
			{
				double[] dp = Problem.InverseProblem(info);
				for (int i = 0; i < dp.Length; i++)
					info.p[i] += dp[i];

				info.V = Problem.DirectProblem(source, receivers, info.p[0]);
			}
		}

		static void Main(string[] args)
		{
			FEMDeltaTest();
		}
	}
}
