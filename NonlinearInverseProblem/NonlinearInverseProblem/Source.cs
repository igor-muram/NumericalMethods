using System;
using System.Numerics;

namespace NonlinearInverseProblem
{
	public class Source
	{
		public Vector3 A { get; set; }
		public Vector3 B { get; set; }
		public double I { get; set; }

		public Source(Vector3 A, Vector3 B, double I)
		{
			this.A = A;
			this.B = B;
			this.I = I;
		}

		public double Potential(Vector3 point, double sigma)
		{
			double rA = Vector3.Distance(A, point);
			double rB = Vector3.Distance(B, point);

			return I * (1.0 / rB - 1.0 / rA) / (2 * Math.PI * sigma);
		}

		public double Potential(Receiver receiver, double sigma) => Potential(receiver.M, sigma) - Potential(receiver.N, sigma);
	}
}
