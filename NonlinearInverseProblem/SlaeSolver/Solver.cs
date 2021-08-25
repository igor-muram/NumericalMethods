using MathUtilities;
using System;

namespace SlaeSolver
{
	public enum SolverTypes { CGM, LOSLU }

	public interface ISolver
	{
		public IMatrix Matrix { get; set; }
		public double[] B { get; set; }
		public int MaxIterCount { get; set; }
		public int IterCount { get; set; }
		public double Eps { get; set; }
		public double Difference { get; set; }

		public double[] Solve();
	}

	public class LOSLU : ISolver
	{
		public IMatrix Matrix { get; set; } = null;
		public double[] B { get; set; } = null;

		public int MaxIterCount { get; set; } = 20000;
		public int IterCount { get; set; } = 0;

		public double Eps { get; set; } = 1.0e-7;
		public double Difference { get; set; } = 0.0;

		double[] Ax { get; set; } = null;
		double[] r { get; set; } = null;
		double[] z { get; set; } = null;
		double[] p { get; set; } = null;
		double[] temp { get; set; } = null;
		int N { get; set; } = 0;

		RawMatrix LU { get; set; }

		struct RawMatrix
		{
			public int N { get; set; }
			public double[] DI { get; set; }
			public double[] AL { get; set; }
			public double[] AU { get; set; }
			public int[] IA { get; set; }
			public int[] JA { get; set; }
		}

		public LOSLU(IMatrix matrix, double[] b)
		{
			N = matrix.N;
			Matrix = matrix;
			B = b;

			Ax = new double[N];
			r = new double[N];
			z = new double[N];
			p = new double[N];
			temp = new double[N];

			LUFactorization();
		}

		public double[] Solve()
		{
			double[] x = new double[N];

			// Calculate r0
			Matrix.Multiply(x, Ax);
			for (int i = 0; i < N; i++)
				r[i] = B[i] - Ax[i];

			Forward(LU, r, r);

			//Calculate z0
			Backward(LU, z, r);

			// Calculate p0
			Matrix.Multiply(z, p);
			Forward(LU, p, p);

			Difference = Utilities.DotProduct(r, r);

			while(IterCount < MaxIterCount && Difference >= Eps)
			{
				// Calculate alpha
				double dotP = Utilities.DotProduct(p, p);
				double a = Utilities.DotProduct(p, r) / dotP;

				// Calculate xk, rk
				for (int i = 0; i < N; i++)
				{
					x[i] += a * z[i];
					r[i] -= a * p[i];
				}

				// Calculate beta
				Backward(LU, Ax, r);
				Matrix.Multiply(Ax, temp);
				Forward(LU, Ax, temp);
				double b = -Utilities.DotProduct(p, Ax) / dotP;

				// Calculate zk, pk
				Backward(LU, temp, r);
				for (int i = 0; i < N; i++)
				{
					z[i] = temp[i] + b * z[i];
					p[i] = Ax[i] + b * p[i];
				}

				// Calculate difference
				Difference = Utilities.DotProduct(r, r);

				IterCount++;
			}

			return x;
		}

		void LUFactorization()
		{
			RawMatrix LU = new RawMatrix();
			LU.N = Matrix.N;
			LU.IA = new int[LU.N + 1];

			for (int i = 0; i < Matrix.N + 1; i++)
				LU.IA[i] = Matrix.IA[i];

			LU.AL = new double[LU.IA[LU.N]];
			LU.AU = new double[LU.IA[LU.N]];
			LU.JA = new int[LU.IA[LU.N]];
			LU.DI = new double[LU.N];

			for (int i = 0; i < Matrix.IA[Matrix.N]; i++)
				LU.JA[i] = Matrix.JA[i];

			for (int i = 0; i < Matrix.N; i++)
			{
				double sumD = 0;
				int i0 = Matrix.IA[i], i1 = Matrix.IA[i + 1];

				for (int k = i0; k < i1; k++)
				{
					double sumL = 0, sumU = 0;
					int j = Matrix.JA[k];

					// Calculate L[i][j], U[j][i]
					int j0 = Matrix.IA[j], j1 = Matrix.IA[j + 1];

					int kl = i0, ku = j0;

					for (; kl < i1 && ku < j1;)
					{
						int j_kl = Matrix.JA[kl];
						int j_ku = Matrix.JA[ku];

						if (j_kl == j_ku)
						{
							sumL += LU.AL[kl] * LU.AU[ku];
							sumU += LU.AU[kl] * LU.AL[ku];
							kl++;
							ku++;
						}
						if (j_kl > j_ku)
							ku++;
						if (j_kl < j_ku)
							kl++;
					}

					LU.AL[k] = Matrix.AL[k] - sumL;
					LU.AU[k] = Matrix.AU[k] - sumU;
					LU.AU[k] /= LU.DI[j];

					// Calculate sum for DI[i]
					sumD += LU.AL[k] * LU.AU[k];
				}

				// Calculate DI[i]
				LU.DI[i] = Matrix.DI[i] - sumD;
			}

			this.LU = LU;
		}

		void Forward(RawMatrix A, double[] x, double[] b)
		{
			double[] di = A.DI;
			double[] al = A.AL;
			double[] au = A.AU;
			int[] ia = A.IA;
			int[] ja = A.JA;
			int N = A.N;


			for (int i = 0; i < N; i++)
			{
				double sum = 0;
				int i0 = ia[i], i1 = ia[i + 1];
				for (int k = i0; k < i1; k++)
				{
					int j = ja[k];
					sum += al[k] * x[j];
				}
				x[i] = (b[i] - sum) / di[i];
			}
		}

		void Backward(RawMatrix A, double[] x, double[] b)
		{
			double[] di = A.DI;
			double[] al = A.AL;
			double[] au = A.AU;
			int[] ia = A.IA;
			int[] ja = A.JA;
			int N = A.N;

			for (int i = 0; i < N; i++)
				x[i] = b[i];

			for (int i = N - 1; i >= 0; i--)
			{
				int i0 = ia[i], i1 = ia[i + 1];
				for (int k = i0; k < i1; k++)
				{
					int j = ja[k];
					x[j] -= au[k] * x[i];

				}
			}
		}
	}

	public class CGM : ISolver
	{
		public IMatrix Matrix { get; set; } = null;
		public double[] B { get; set; } = null;

		public int MaxIterCount { get; set; } = 200000;
		public int IterCount { get; set; } = 0;

		public double Eps { get; set; } = 1.0e-15;
		public double Difference { get; set; } = 0.0;

		double[] Ax { get; set; } = null;
		double[] r { get; set; } = null;
		double[] z { get; set; } = null;
		int N { get; set; } = 0;

		public CGM(IMatrix matrix, double[] b)
		{
			N = matrix.N;
			Matrix = matrix;
			B = b;

			Ax = new double[N];
			r = new double[N];
			z = new double[N];
		}

		public double[] Solve()
		{
			double[] x = new double[N];

			// Calculate r0, z0
			Matrix.Multiply(x, Ax);
			for (int i = 0; i < N; i++)
			{
				r[i] = B[i] - Ax[i];
				z[i] = r[i];
			}

			double dotF = Utilities.DotProduct(B, B);
			double dotR0 = Utilities.DotProduct(r, r);
			Difference = dotR0 / dotF;

			while(IterCount < MaxIterCount && Difference >= Eps * Eps)
			{
				// Calculate alpha
				Matrix.Multiply(z, Ax);
				double a = dotR0 / Utilities.DotProduct(Ax, z);

				// Calculate xk, rk
				for (int i = 0; i < N; i++)
				{
					x[i] += a * z[i];
					r[i] -= a * Ax[i];
				}

				// Calculate beta
				double dotR1 = Utilities.DotProduct(r, r);
				double b = dotR1 / dotR0;

				// Calculate zk
				for (int i = 0; i < N; i++)
					z[i] = r[i] + b * z[i];

				// Calculate difference
				Difference = dotR1 / dotF;
				dotR0 = dotR1;
			}

			return x;
		}
	}
}
