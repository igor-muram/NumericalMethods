using System;
using System.Linq;

namespace SlaeSolver
{

	public interface IMatrix
	{ 
		public int N { get; set; }
		public int[] IA { get; set; }
		public int[] JA { get; set; }
		public double[] DI { get; set; }
		public double[] AL { get; set; }
		public double[] AU { get; set; }

		public void Set(int i, int j, double value);
		public void Multiply(double[] vector, double[] result);
	}

	public class FullSparseMatrix : IMatrix
	{
		public int N { get; set; }
		public int[] IA { get; set; }
		public int[] JA { get; set; }
		public double[] DI { get; set; }
		public double[] AL { get; set; }
		public double[] AU { get; set; }

		public FullSparseMatrix(int N)
		{
			this.N = N;
			DI = new double[N];
			IA = new int[N + 1];
			JA = new int[N * (N - 1) / 2];
			AL = new double[N * (N - 1) / 2];
			AU = new double[N * (N - 1) / 2];

			for (int i = 0; i < N; i++)
				IA[i + 1] = IA[i] + i;

			int k = 0;
			for (int i = 1; i < N; i++)
				for (int j = 0; j < i; j++)
				{
					JA[k] = j;
					k++;
				}
		}

		public void Set(int i, int j, double value)
		{
			if (i >= N || j >= N || i < 0 || j < 0)
				throw new ArgumentOutOfRangeException("Bad index for matrix");

			if (i == j)
				DI[i] = value;

			if (j > i)
			{
				int width = IA[j + 1] - IA[j];
				int begin = j - width;
				int end = j;

				for (int k = begin; k < end; k++)
					if (k == i)
						AU[IA[j] + k - begin] = value;
			}
			else
			{
				int width = IA[i + 1] - IA[i];
				int begin = i - width;
				int end = i;

				for (int k = begin; k < end; k++)
					if (k == j)
						AL[IA[i] + k - begin] = value;
			}
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
					result[j] += AU[k] * vector[i];
				}
			}
		}
	}

	public class SparseMatrix : IMatrix
	{
		public int N { get; set; }
		public int[] IA { get; set; }
		public int[] JA { get; set; }
		public double[] DI { get; set; }
		public double[] AL { get; set; }
		public double[] AU { get; set; }

		public SparseMatrix(int N, MatrixPortrait portrait)
		{
			this.N = N;
			DI = new double[N];
			IA = portrait.IA;
			JA = portrait.JA;
			AL = new double[N * (N - 1) / 2];
			AU = new double[N * (N - 1) / 2];
		}

		public void Set(int i, int j, double value)
		{
			if (i >= N || j >= N || i < 0 || j < 0)
				throw new ArgumentOutOfRangeException("Bad index for matrix");

			if (i == j)
				DI[i] = value;

			if (j > i)
			{
				int width = IA[j + 1] - IA[j];
				int begin = j - width;
				int end = j;

				for (int k = begin; k < end; k++)
					if (k == i)
						AU[IA[j] + k - begin] = value;
			}
			else
			{
				int width = IA[i + 1] - IA[i];
				int begin = i - width;
				int end = i;

				for (int k = begin; k < end; k++)
					if (k == j)
						AL[IA[i] + k - begin] = value;
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
					result[j] += AU[k] * vector[i];
				}
			}
		}
	}

	public class SymmetricSparseMatrix : IMatrix
	{
		public int N { get; set; }
		public int[] IA { get; set; }
		public int[] JA { get; set; }
		public double[] DI { get; set; }
		public double[] AL { get; set; }
		public double[] AU { get => AL; set => throw new NotImplementedException(); }

		public SymmetricSparseMatrix(int N, MatrixPortrait portrait)
		{
			this.N = N;
			DI = new double[N];
			IA = portrait.IA;
			JA = portrait.JA;
			AL = new double[portrait.IA.Last()];
		}

		public void Set(int i, int j, double value)
		{
			if (i >= N || j >= N || i < 0 || j < 0)
				throw new ArgumentOutOfRangeException("Bad index for matrix");

			if (i == j)
				DI[i] = value;

			if (j > i)
			{
				int width = IA[j + 1] - IA[j];
				int begin = j - width;
				int end = j;

				for (int k = begin; k < end; k++)
					if (k == i)
						AL[IA[j] + k - begin] = value;
			}
			else
			{
				int width = IA[i + 1] - IA[i];
				int begin = i - width;
				int end = i;

				for (int k = begin; k < end; k++)
					if (k == j)
						AL[IA[i] + k - begin] = value;
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
	}
}
