using System.Collections;
using System.Collections.Generic;

namespace FEM
{
	public class Mesh : IEnumerable
	{
		public List<FiniteElement> Elements { get; set; } = new List<FiniteElement>();

		public IEnumerator GetEnumerator()
		{
			return Elements.GetEnumerator();
		}
	}

	public class MeshBuilder
	{
		public Dictionary<int, int>[] EdgeMatrix { get; set; }
		public int NodeCount { get; set; } = 0;

		public MeshBuilder(int nodeCount)
		{
			EdgeMatrix = new Dictionary<int, int>[nodeCount];
			for (int i = 0; i < nodeCount; i++)
				EdgeMatrix[i] = new Dictionary<int, int>();

			NodeCount = nodeCount;
		}

		public int Build(Mesh mesh)
		{
			foreach (FiniteElement e in mesh)
			{
				int index = 3;
				for (int i = 0; i < 3; i++)
					for (int j = i + 1; j < 3; j++, index += 2)
					{
						int a = e.Vertices[i];
						int b = e.Vertices[j];
						bool f = a > b;
						if (f) (a, b) = (b, a);

						if (!EdgeMatrix[a].ContainsKey(b))
						{
							e.Vertices[index] = NodeCount + (f ? 1 : 0);
							e.Vertices[index + 1] = NodeCount + (f ? 0 : 1);
							EdgeMatrix[a][b] = NodeCount;
							NodeCount += 2;
						}
						else
						{
							int a1 = EdgeMatrix[a][b];
							e.Vertices[index] = a1 + (f ? 1 : 0);
							e.Vertices[index + 1] = a1 + (f ? 0 : 1);
						}
					}
				e.Vertices[index] = NodeCount;
				NodeCount++;
			}

			return NodeCount;
		}

		public void BuildBoundary(FirstBoundary boundary)
		{
			foreach (Edge edge in boundary.Edges)
			{
				int a = edge.V1;
				int b = edge.V4;
				bool f = a > b;
				if (f) (a, b) = (b, a);

				edge.V2 = EdgeMatrix[a][b] + (f ? 1 : 0);
				edge.V3 = EdgeMatrix[a][b] + (f ? 0 : 1);
			}
		}
	}
}
