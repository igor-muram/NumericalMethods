using System.Numerics;

namespace NonlinearInverseProblem
{
	public class Receiver
	{
		public Vector3 M { get; set; }
		public Vector3 N { get; set; }

		public Receiver(Vector3 M, Vector3 N)
		{
			this.M = M;
			this.N = N;
		}
	}
}
