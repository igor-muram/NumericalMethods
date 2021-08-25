using MathUtilities;
using SlaeSolver;
using System;
using System.Collections.Generic;
using System.Reflection.Metadata.Ecma335;

namespace FEM
{
	public struct ProblemInfo
	{
		public Point[] Points { get; set; }
		public Dictionary<int, Material> Materials { get; set; }
		public Mesh Mesh { get; set; }
		public FirstBoundary FB { get; set; }

		public Func<double, double, double> F { get; set; }

		public SolverTypes SolverType { get; set; }
	}

	public class FEMrz
	{
		public ProblemInfo Info { get; set; }
		public ISolver Solver { get; set; } = null;

		int NodeCount { get; set; } = 0;
		MeshBuilder mb { get; set; } = null;
		PortraitBuilder pb { get; set; } = null;
		SLAEBuilder sb { get; set; } = null;

		public double[] Weights { get; set; } = null;

		public FEMrz(ProblemInfo info)
		{
			Info = info;

			switch (Info.SolverType)
			{
				case SolverTypes.LOSLU:
					Solver = new LOSLU();
					break;

				default:
					Solver = new LOSLU();
					break;
			}
		}

		public void Solve()
		{
			mb = new MeshBuilder(Info.Points.Length);
			NodeCount = mb.Build(Info.Mesh);

			mb.BuildBoundary(Info.FB);

			pb = new PortraitBuilder(NodeCount, Info.Mesh);
			MatrixPortrait mp = new MatrixPortrait();
			pb.Build(mp);

			IMatrix A = new SymmetricSparseMatrix(NodeCount, mp);
			double[] B = new double[NodeCount];

			SLAEInfo slaeinfo = new SLAEInfo();
			slaeinfo.Points = Info.Points;
			slaeinfo.Mesh = Info.Mesh;
			slaeinfo.F = Info.F;

			sb = new SLAEBuilder(slaeinfo);
			sb.Build(A, B);
			sb.AddBoundary(A, B, Info.FB);

			Weights = Solver.Solve(A, B);
		}

		public double U(Point p)
		{
			FiniteElement element = null;

			foreach (FiniteElement e in Info.Mesh)
			{
				if (Utilities.PointInsideTriangle(Info.Points[e.V1], Info.Points[e.V2], Info.Points[e.V3], p))
				{
					element = e;
					break;
				}
			}

			if (element != null)
			{
				(double L1, double L2, double L3) = Utilities.GetL(Info.Points[element.V1], Info.Points[element.V2], Info.Points[element.V3], p);

				double result = 0.0;
				for (int i = 0; i < FEMInfo.BasisSize; i++)
					result += element.Vertices[i] * FEMInfo.LBasis[i](L1, L2, L3);

				return result;				
			}

			return 0.0;
		}
	}

	public struct ProblemInfoDelta
	{
		public Point[] Points { get; set; }
		public Dictionary<int, Material> Materials { get; set; }
		public Mesh Mesh { get; set; }
		public FirstBoundary FB { get; set; }

		public double R0 { get; set; }
		public double Z0 { get; set; }
		public double DeltaCoefficient { get; set; }

		public SolverTypes SolverType { get; set; }
	}

	public class FEMrzDelta
	{
		public ProblemInfoDelta Info { get; set; }
		public ISolver Solver { get; set; } = null;

		int NodeCount { get; set; } = 0;
		MeshBuilder mb { get; set; } = null;
		PortraitBuilder pb { get; set; } = null;
		SLAEBuilderDelta sb { get; set; } = null;

		public double[] Weights { get; set; } = null;

		public FEMrzDelta(ProblemInfoDelta info)
		{
			Info = info;
			switch (Info.SolverType)
			{
				case SolverTypes.LOSLU:
					Solver = new LOSLU();
					break;

				default:
					Solver = new LOSLU();
					break;
			}
		}

		public void Solve()
		{
			mb = new MeshBuilder(Info.Points.Length);
			NodeCount = mb.Build(Info.Mesh);

			mb.BuildBoundary(Info.FB);

			pb = new PortraitBuilder(NodeCount, Info.Mesh);
			MatrixPortrait mp = new MatrixPortrait();
			pb.Build(mp);

			IMatrix A = new SymmetricSparseMatrix(NodeCount, mp);
			double[] B = new double[NodeCount];

			SLAEInfoDelta slaeinfo = new SLAEInfoDelta();
			slaeinfo.Points = Info.Points;
			slaeinfo.Mesh = Info.Mesh;
			slaeinfo.R0 = Info.R0;
			slaeinfo.Z0 = Info.Z0;
			slaeinfo.DeltaCoefficient = Info.DeltaCoefficient;

			sb = new SLAEBuilderDelta(slaeinfo);
			sb.Build(A, B);
			sb.AddBoundary(A, B, Info.FB);

			Console.WriteLine("Start Solving SLAE");
			Weights = Solver.Solve(A, B);
		}

		public double U(Point p)
		{
			FiniteElement element = null;

			foreach (FiniteElement e in Info.Mesh)
			{
				if (Utilities.PointInsideTriangle(Info.Points[e.V1], Info.Points[e.V2], Info.Points[e.V3], p))
				{
					element = e;
					break;
				}
			}

			if (element != null)
			{
				(double L1, double L2, double L3) = Utilities.GetL(Info.Points[element.V1], Info.Points[element.V2], Info.Points[element.V3], p);

				double result = 0.0;
				for (int i = 0; i < FEMInfo.BasisSize; i++)
					result += Weights[element.Vertices[i]] * FEMInfo.LBasis[i](L1, L2, L3);

				return result;
			}

			return 0.0;
		}
	}
}
