using System.Collections.Generic;

namespace FEM
{
	public class FiniteElement
	{
		public int[] Vertices { get; set; } = new int[FEMInfo.BasisSize];
		public Material Material { get; set; } = null;

		public int V1 => Vertices[0];
		public int V2 => Vertices[1];
		public int V3 => Vertices[2];
	}
}
