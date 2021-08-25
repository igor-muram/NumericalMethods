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
	}
}
