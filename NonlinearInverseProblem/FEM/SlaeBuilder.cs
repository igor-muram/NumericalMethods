using MathUtilities;
using SlaeSolver;
using System;
using System.Collections.Generic;

namespace FEM
{
	public struct SLAEInfo
	{
		public Point[] Points { get; set; }
		public Mesh Mesh { get; set; }
		public Func<double, double, double> F { get; set; }
	}

	public class SLAEBuilder
	{
		SLAEInfo Info;

		double[,] local;
		double[] localb;

		List<FEMInfo.LocalCompG>[,] gPattern;
		List<FEMInfo.LocalCompM>[,] mPattern;

		double[][] M;

		public SLAEBuilder(SLAEInfo info)
		{
			Info = info;

			local = new double[FEMInfo.BasisSize, FEMInfo.BasisSize];
			localb = new double[FEMInfo.BasisSize];

			BuildPatterns();
		}

		public void Build(IMatrix A, double[] b)
		{
			foreach (FiniteElement e in Info.Mesh)
			{
				BuildLocalMatrix(e);
				BuildLocalB(e);
				AddLocalToGlobal(A, b, e);
			}
		}

		public void AddBoundary(IMatrix A, double[] b, FirstBoundary boundary)
		{
			Func<double, double>[] basis = new Func<double, double>[4]
			{
				(double x) => 0.5 * (1 - x) * (3 * (1 - x) - 1) * (3 * (1 - x) - 2),
				(double x) => 4.5 * (1 - x) * x * (3 * (1 - x) - 1),
				(double x) => 4.5 * (1 - x) * x * (3 * x - 1),
				(double x) => 0.5 * x * (3 * x - 1) * (3 * x - 2)
			};

			foreach (Edge edge in boundary.Edges)
			{
				double[][] M = new double[4][]
				{
					new double[4] { 8.0 / 105,    33.0 / 560,    -3.0 / 140,    19.0 / 1680 },
					new double[4] { 33.0 / 560,   27.0 / 70,     -27.0 / 560,   -3.0 / 140 },
					new double[4] { -3.0 / 140,   -27.0 / 560,   27.0 / 70,     33.0 / 560 },
					new double[4] { 19.0 / 1680,  -3.0 / 140,    33.0 / 560,    8.0 / 105 }
				};

				double[] f = new double[4];

				double r0 = Info.Points[edge.V1].R;
				double r1 = Info.Points[edge.V4].R;
				double z0 = Info.Points[edge.V1].Z;
				double z1 = Info.Points[edge.V4].Z;

				for (int i = 0; i < 4; i++)
					f[i] = Utilities.NewtonCotes(0.0, 1.0, (double t) => edge.Function(r0 + t * (r1 - r0), z0 + t * (z1 - z0)) * basis[i](t));

				double[] q = Utilities.Gauss(M, f);

				A.DI[edge.V1] = 1.0E+50;
				b[edge.V1] = 1.0E+50 * q[0];
				A.DI[edge.V2] = 1.0E+50;
				b[edge.V2] = 1.0E+50 * q[1];
				A.DI[edge.V3] = 1.0E+50;
				b[edge.V3] = 1.0E+50 * q[2];
				A.DI[edge.V4] = 1.0E+50;
				b[edge.V4] = 1.0E+50 * q[3];
			}
		}

		void BuildPatterns()
		{
			gPattern = new List<FEMInfo.LocalCompG>[FEMInfo.BasisSize, FEMInfo.BasisSize];
			mPattern = new List<FEMInfo.LocalCompM>[FEMInfo.BasisSize, FEMInfo.BasisSize];

			for (int i = 0; i < FEMInfo.BasisSize; i++)
				for (int j = 0; j < FEMInfo.BasisSize; j++)
				{
					gPattern[i, j] = new List<FEMInfo.LocalCompG>();
					mPattern[i, j] = new List<FEMInfo.LocalCompM>();
				}

			M = new double[FEMInfo.BasisSize][];
			for (int i = 0; i < FEMInfo.BasisSize; i++)
				M[i] = new double[FEMInfo.BasisSize];

			for (int k = 0; k < 3; k++)
			{
				int rnum = 0;
				int v11 = 0, v22 = 0, v33 = 0;
				if (k == 0) { v11++; rnum = 0; }
				else if (k == 1) { v22++; rnum = 1; }
				else { v33++; rnum = 2; }

				for (int i = 0; i < FEMInfo.BasisSize; i++)
					for (int j = 0; j < FEMInfo.BasisSize; j++)
					{
						foreach (var r in FEMInfo.Ders[i])
							foreach (var z in FEMInfo.Ders[j])
							{
								int v1 = r.V1 + z.V1 + v11;
								int v2 = r.V2 + z.V2 + v22;
								int v3 = r.V3 + z.V3 + v33;

								double coeff = r.Coefficient * z.Coefficient * Utilities.Factorial(v1) * Utilities.Factorial(v2) * Utilities.Factorial(v3) / (double)Utilities.Factorial(v1 + v2 + v3 + 2);
								gPattern[i, j].Add(new FEMInfo.LocalCompG(r.GradNo, z.GradNo, coeff, rnum));
							}

						double sum = 0;
						foreach (var r in FEMInfo.Basis[i])
							foreach (var z in FEMInfo.Basis[j])
							{
								int v1 = r.V1 + z.V1 + v11;
								int v2 = r.V2 + z.V2 + v22;
								int v3 = r.V3 + z.V3 + v33;

								sum += r.Coefficient * z.Coefficient * Utilities.Factorial(v1) * Utilities.Factorial(v2) * Utilities.Factorial(v3) / (double)Utilities.Factorial(v1 + v2 + v3 + 2);
							}

						mPattern[i, j].Add(new FEMInfo.LocalCompM(sum, rnum));
					}
			}
		}

		void BuildLocalMatrix(FiniteElement e)
		{
			Point a = Info.Points[e.V1];
			Point b = Info.Points[e.V2];
			Point c = Info.Points[e.V3];

			double D = Math.Abs(Utilities.Det(a, b, c));
			double[] alpha = Utilities.Alpha(a, b, c);

			Grad[] grads = new Grad[3]
			{
				new Grad(alpha[0], alpha[1]),
				new Grad(alpha[2], alpha[3]),
				new Grad(alpha[4], alpha[5])
			};

			for (int i = 0; i < FEMInfo.BasisSize; i++)
				Array.Fill(M[i], 0.0);

			for (int i = 0; i < FEMInfo.BasisSize; i++)
				for (int j = 0; j < FEMInfo.BasisSize; j++)
				{
					double G = 0;
					foreach (FEMInfo.LocalCompG comp in gPattern[i, j])
					{
						double scalGrad = grads[comp.Grad1].R * grads[comp.Grad2].R + grads[comp.Grad1].Z * grads[comp.Grad2].Z;
						G += comp.Coefficient * scalGrad * Info.Points[e.Vertices[comp.Rnum]].R;
					}

					foreach (FEMInfo.LocalCompM comp in mPattern[i, j])
						M[i][j] += comp.Coefficient * Info.Points[e.Vertices[comp.Rnum]].R;

					local[i, j] = (e.Material.Gamma * M[i][j] + e.Material.Sigma * G) * D;
				}
		}

		void BuildLocalB(FiniteElement e)
		{
			Point a = Info.Points[e.V1];
			Point b = Info.Points[e.V2];
			Point c = Info.Points[e.V3];

			Point[] coords = Utilities.TriangleCubicPoints(a, b, c);
			double D = Math.Abs(Utilities.Det(a, b, c));

			double[] temp = new double[FEMInfo.BasisSize];
			for (int i = 0; i < FEMInfo.BasisSize; i++)
				temp[i] = Info.F(coords[i].R, coords[i].Z);

			for (int i = 0; i < FEMInfo.BasisSize; i++)
			{
				localb[i] = temp[0] * M[i][0];

				for (int j = 1; j < FEMInfo.BasisSize; j++)
					localb[i] += temp[j] * M[i][j];

				localb[i] *= D;
			}
		}

		void AddLocalToGlobal(IMatrix A, double[] B, FiniteElement e)
		{
			for (int i = 0; i < FEMInfo.BasisSize; i++)
			{
				A.DI[e.Vertices[i]] += local[i, i];
				B[e.Vertices[i]] += localb[i];

				for (int j = 0; j < i; j++)
				{
					int a = e.Vertices[i];
					int b = e.Vertices[j];
					if (a < b) (a, b) = (b, a);

					if (A.IA[a + 1] > A.IA[a])
					{
						Span<int> span = new Span<int>(A.JA, A.IA[a], A.IA[a + 1] - A.IA[a]);

						int k;
						for (k = 0; k < A.IA[a + 1] - A.IA[a]; k++)
							if (span[k] == b)
								break;

						int index = A.IA[a] + k;
						A.AL[index] += local[i, j];
					}
				}
			}
		}
	};

	public struct SLAEInfoDelta
	{
		public Point[] Points { get; set; }
		public Mesh Mesh { get; set; }

		public double R0 { get; set; }
		public double Z0 { get; set; }

		public double DeltaCoefficient { get; set; }
	}

	public class SLAEBuilderDelta
	{
		SLAEInfoDelta Info;

		double[,] local;
		double[] localb;

		List<FEMInfo.LocalCompG>[,] gPattern;
		List<FEMInfo.LocalCompM>[,] mPattern;

		double[][] M;

		public SLAEBuilderDelta(SLAEInfoDelta info)
		{
			Info = info;

			local = new double[FEMInfo.BasisSize, FEMInfo.BasisSize];
			localb = new double[FEMInfo.BasisSize];

			BuildPatterns();
		}

		public void Build(IMatrix A, double[] b)
		{
			foreach (FiniteElement e in Info.Mesh)
			{
				BuildLocalMatrix(e);
				BuildLocalB(e);
				AddLocalToGlobal(A, b, e);
			}
		}

		public void AddBoundary(IMatrix A, double[] b, FirstBoundary boundary)
		{
			Func<double, double>[] basis = new Func<double, double>[4]
			{
				(double x) => 0.5 * (1 - x) * (3 * (1 - x) - 1) * (3 * (1 - x) - 2),
				(double x) => 4.5 * (1 - x) * x * (3 * (1 - x) - 1),
				(double x) => 4.5 * (1 - x) * x * (3 * x - 1),
				(double x) => 0.5 * x * (3 * x - 1) * (3 * x - 2)
			};

			foreach (Edge edge in boundary.Edges)
			{
				double[][] M = new double[4][]
				{
					new double[4] { 8.0 / 105,    33.0 / 560,    -3.0 / 140,    19.0 / 1680 },
					new double[4] { 33.0 / 560,   27.0 / 70,     -27.0 / 560,   -3.0 / 140 },
					new double[4] { -3.0 / 140,   -27.0 / 560,   27.0 / 70,     33.0 / 560 },
					new double[4] { 19.0 / 1680,  -3.0 / 140,    33.0 / 560,    8.0 / 105 }
				};

				double[] f = new double[4];

				double r0 = Info.Points[edge.V1].R;
				double r1 = Info.Points[edge.V4].R;
				double z0 = Info.Points[edge.V1].Z;
				double z1 = Info.Points[edge.V4].Z;

				for (int i = 0; i < 4; i++)
					f[i] = Utilities.NewtonCotes(0.0, 1.0, (double t) => edge.Function(r0 + t * (r1 - r0), z0 + t * (z1 - z0)) * basis[i](t));

				double[] q = Utilities.Gauss(M, f);

				if (!(Math.Abs(r0 - Info.R0) < 1.0e-15 && Math.Abs(z0 - Info.Z0) < 1.0e-15))
				{
					A.DI[edge.V1] = 1.0E+50;
					b[edge.V1] = 1.0E+50 * q[0];
				}
				A.DI[edge.V2] = 1.0E+50;
				b[edge.V2] = 1.0E+50 * q[1];
				A.DI[edge.V3] = 1.0E+50;
				b[edge.V3] = 1.0E+50 * q[2];

				if (!(Math.Abs(r1 - Info.R0) < 1.0e-15 && Math.Abs(z1 - Info.Z0) < 1.0e-15))
				{
					A.DI[edge.V4] = 1.0E+50;
					b[edge.V4] = 1.0E+50 * q[3];
				}
			}
		}

		void BuildPatterns()
		{
			gPattern = new List<FEMInfo.LocalCompG>[FEMInfo.BasisSize, FEMInfo.BasisSize];
			mPattern = new List<FEMInfo.LocalCompM>[FEMInfo.BasisSize, FEMInfo.BasisSize];

			for (int i = 0; i < FEMInfo.BasisSize; i++)
				for (int j = 0; j < FEMInfo.BasisSize; j++)
				{
					gPattern[i, j] = new List<FEMInfo.LocalCompG>();
					mPattern[i, j] = new List<FEMInfo.LocalCompM>();
				}

			M = new double[FEMInfo.BasisSize][];
			for (int i = 0; i < FEMInfo.BasisSize; i++)
				M[i] = new double[FEMInfo.BasisSize];

			for (int k = 0; k < 3; k++)
			{
				int rnum = 0;
				int v11 = 0, v22 = 0, v33 = 0;
				if (k == 0) { v11++; rnum = 0; }
				else if (k == 1) { v22++; rnum = 1; }
				else { v33++; rnum = 2; }

				for (int i = 0; i < FEMInfo.BasisSize; i++)
					for (int j = 0; j < FEMInfo.BasisSize; j++)
					{
						foreach (var r in FEMInfo.Ders[i])
							foreach (var z in FEMInfo.Ders[j])
							{
								int v1 = r.V1 + z.V1 + v11;
								int v2 = r.V2 + z.V2 + v22;
								int v3 = r.V3 + z.V3 + v33;

								double coeff = r.Coefficient * z.Coefficient * Utilities.Factorial(v1) * Utilities.Factorial(v2) * Utilities.Factorial(v3) / (double)Utilities.Factorial(v1 + v2 + v3 + 2);
								gPattern[i, j].Add(new FEMInfo.LocalCompG(r.GradNo, z.GradNo, coeff, rnum));
							}

						double sum = 0;
						foreach (var r in FEMInfo.Basis[i])
							foreach (var z in FEMInfo.Basis[j])
							{
								int v1 = r.V1 + z.V1 + v11;
								int v2 = r.V2 + z.V2 + v22;
								int v3 = r.V3 + z.V3 + v33;

								sum += r.Coefficient * z.Coefficient * Utilities.Factorial(v1) * Utilities.Factorial(v2) * Utilities.Factorial(v3) / (double)Utilities.Factorial(v1 + v2 + v3 + 2);
							}

						mPattern[i, j].Add(new FEMInfo.LocalCompM(sum, rnum));
					}
			}
		}

		void BuildLocalMatrix(FiniteElement e)
		{
			Point a = Info.Points[e.V1];
			Point b = Info.Points[e.V2];
			Point c = Info.Points[e.V3];

			double D = Math.Abs(Utilities.Det(a, b, c));
			double[] alpha = Utilities.Alpha(a, b, c);

			Grad[] grads = new Grad[3]
			{
				new Grad(alpha[0], alpha[1]),
				new Grad(alpha[2], alpha[3]),
				new Grad(alpha[4], alpha[5])
			};

			for (int i = 0; i < FEMInfo.BasisSize; i++)
				Array.Fill(M[i], 0.0);

			for (int i = 0; i < FEMInfo.BasisSize; i++)
				for (int j = 0; j < FEMInfo.BasisSize; j++)
				{
					double G = 0;
					foreach (FEMInfo.LocalCompG comp in gPattern[i, j])
					{
						double scalGrad = grads[comp.Grad1].R * grads[comp.Grad2].R + grads[comp.Grad1].Z * grads[comp.Grad2].Z;
						G += comp.Coefficient * scalGrad * Info.Points[e.Vertices[comp.Rnum]].R;
					}

					foreach (FEMInfo.LocalCompM comp in mPattern[i, j])
						M[i][j] += comp.Coefficient * Info.Points[e.Vertices[comp.Rnum]].R;

					local[i, j] = (e.Material.Gamma * M[i][j] + e.Material.Sigma * G) * D;
				}
		}

		void BuildLocalB(FiniteElement e)
		{
			Point a = Info.Points[e.V1];
			Point b = Info.Points[e.V2];
			Point c = Info.Points[e.V3];
			Point[] coords = Utilities.TriangleCubicPoints(a, b, c);

			for (int i = 0; i < FEMInfo.BasisSize; i++)
			{
				if (Math.Abs(coords[i].R - Info.R0) < 1.0e-15 && Math.Abs(coords[i].Z - Info.Z0) < 1.0e-15)
					localb[i] = 1.0 * Info.DeltaCoefficient;
				else
					localb[i] = 0.0;
			}
		}

		void AddLocalToGlobal(IMatrix A, double[] B, FiniteElement e)
		{
			for (int i = 0; i < FEMInfo.BasisSize; i++)
			{
				A.DI[e.Vertices[i]] += local[i, i];
				B[e.Vertices[i]] += localb[i];

				for (int j = 0; j < i; j++)
				{
					int a = e.Vertices[i];
					int b = e.Vertices[j];
					if (a < b) (a, b) = (b, a);

					if (A.IA[a + 1] > A.IA[a])
					{
						Span<int> span = new Span<int>(A.JA, A.IA[a], A.IA[a + 1] - A.IA[a]);

						int k;
						for (k = 0; k < A.IA[a + 1] - A.IA[a]; k++)
							if (span[k] == b)
								break;

						int index = A.IA[a] + k;
						A.AL[index] += local[i, j];
					}
				}
			}
		}
	};
}