using System;
using System.Linq;

namespace ElectromagneticProblem.Matrix
{
	public interface IMatrix
	{
		public int N { get; set; }
		MatrixPortrait Portrait { get; set; }
		public double[] DI { get; set; }
		public double[] AL { get; set; }
		public double[] AU { get; set; }
		public int[] IA => Portrait.IA;
		public int[] JA => Portrait.JA;


		public void Set(int i, int j, double value);
		public void Add(int i, int j, double value);
		public void Multiply(double[] vector, double[] result);
		public void Clear();
	}

	public class SymmetricSparseMatrix : IMatrix
	{
		public int N { get; set; }
		public MatrixPortrait Portrait { get; set; }

		public int[] IA => Portrait.IA;
		public int[] JA => Portrait.JA;
		public double[] DI { get; set; }
		public double[] AL { get; set; }
		public double[] AU { get => AL; set => throw new NotImplementedException(); }

		public SymmetricSparseMatrix(int N, MatrixPortrait portrait)
		{
			this.N = N;
			Portrait = portrait;

			DI = new double[N];
			AL = new double[portrait.IA.Last()];
		}

		public void Set(int i, int j, double value)
		{
			if (i >= N || j >= N || i < 0 || j < 0)
				throw new ArgumentOutOfRangeException("Bad index for matrix");

			if (i == j)
			{
				DI[i] = value;
				return;
			}

			if (j > i)
				(i, j) = (j, i);

			for (int k = IA[i]; k < IA[i + 1]; k++)
				if (JA[k] == j)
				{
					AL[k] = value;
					return;
				}

			throw new ArgumentOutOfRangeException("Bad index for matrix");
		}

		public void Add(int i, int j, double value)
		{
			if (i >= N || j >= N || i < 0 || j < 0)
				throw new ArgumentOutOfRangeException("Bad index for matrix");

			if (i == j)
			{
				DI[i] += value;
				return;
			}

			if (j > i)
				(i, j) = (j, i);

			for (int k = IA[i]; k < IA[i + 1]; k++)
				if (JA[k] == j)
				{
					AL[k] += value;
					return;
				}

			throw new ArgumentOutOfRangeException("Bad index for matrix");
		}

		public void Multiply(double[] vector, double[] result)
		{
			if (vector.Length != N)
				throw new Exception("Vector's length is not equal Matrix's size");

			for (int i = 0; i < N; i++)
			{
				result[i] = vector[i] * DI[i];
				for (int k = IA[i]; k < IA[i + 1]; k++)
				{
					int j = JA[k];
					result[i] += AL[k] * vector[j];
					result[j] += AL[k] * vector[i];
				}
			}
		}

		public void Clear()
		{
			Array.Fill(DI, 0.0);
			Array.Fill(AL, 0.0);
			Array.Fill(AU, 0.0);
		}
	}
}
