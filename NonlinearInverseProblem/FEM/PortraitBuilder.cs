using System.Collections.Generic;
using SlaeSolver;

namespace FEM
{
	public class PortraitBuilder
	{
		public int NodeCount { get; set; }
		public SortedSet<int>[] connections { get; set; }

		public PortraitBuilder(int nodeCount, Mesh mesh)
		{
			NodeCount = nodeCount;
			connections = new SortedSet<int>[nodeCount];
			for (int i = 0; i < NodeCount; i++)
				connections[i] = new SortedSet<int>();

			BuildConnections(mesh);
		}

		public void Build(MatrixPortrait MP)
		{
			MP.IA = new int[NodeCount + 1];
			MP.IA[0] = MP.IA[1] = 0;

			for (int i = 2; i <= NodeCount; i++)
			{
				int col = MP.IA[i - 1];
				MP.IA[i] = col + connections[i - 1].Count;
			}

			MP.JA = new int[MP.IA[NodeCount]];

			for (int i = 1, k = 0; i < NodeCount; i++)
			{
				foreach (int j in connections[i])
				{
					MP.JA[k] = j;
					k++;
				}
			}
		}

		void BuildConnections(Mesh mesh)
		{
			foreach (FiniteElement e in mesh)
				for (int i = 1; i < FEMInfo.BasisSize; i++)
					for (int j = 0; j < i; j++)
					{
						int a = e.Vertices[i];
						int b = e.Vertices[j];
						if (a < b) (a, b) = (b, a);

						connections[a].Add(b);
					}
		}
	};
}
