using System;
using System.Collections.Generic;

namespace FEM
{
	public enum ConditionType { First, SecondNull }

	public class Edge
	{
		public int V1 { get; set; } = 0;
		public int V2 { get; set; } = 0;
		public int V3 { get; set; } = 0;
		public int V4 { get; set; } = 0;
		public Func<double, double, double> Function { get; set; } = null;
	}

	public class FirstBoundary
	{
		public List<Edge> Edges { get; set; } = new List<Edge>();
	}
}
