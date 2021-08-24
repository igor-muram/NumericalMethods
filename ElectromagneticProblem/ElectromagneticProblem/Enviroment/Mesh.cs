using MathUtility;
using System.Collections.Generic;

namespace ElectromagneticProblem.Enviroment
{
	public class FiniteElement
	{
		public int[] Vertices;
		public Material Material;

		public FiniteElement(int verticesNum)
		{
			Vertices = new int[verticesNum];
			Material = new Material();
		}

		public FiniteElement(int verticesNum, Material material)
		{
			Vertices = new int[verticesNum];
			Material = material;
		}

		public int this[int i]
		{
			get { return Vertices[i]; }
			set { Vertices[i] = value; }
		}
	}

	public class Mesh
	{
		public int NodeCount { get; set; } = 0;
		public Point[] Points { get; set; } = null;
		public List<FiniteElement> Elements { get; set; } = new List<FiniteElement>();
		public FirstNullBoundary FirstBoundary { get; set; } = null;
	}
}
