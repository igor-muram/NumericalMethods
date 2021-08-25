using MathUtilities;
using System;
using System.Collections.Generic;

namespace FEM
{
	public enum AreaSide { Left, Top, Right, Bottom }

	public struct AreaInfo
	{
		public double R0 { get; set; }
		public double Z0 { get; set; }
		public double Width { get; set; }

		public double FirstLayerHeight { get; set; }
		public double SecondLayerHeight { get; set; }

		public int HorizontalNodeCount { get; set; }
		public int FirstLayerVerticalNodeCount { get; set; }
		public int SecondLayerVerticalNodeCount { get; set; }

		public Dictionary<int, Material> Materials { get; set; }
		public Dictionary<AreaSide, (ConditionType, Func<double, double, double>)> Conditions { get; set; }

	}

	public class GridBuilder
	{
		public AreaInfo Info { get; set; }

		public Mesh Grid { get; set; } = null;
		public List<Point> Points { get; set; } = null;
		public FirstBoundary FB { get; set; } = null;

		public GridBuilder(AreaInfo info)
		{
			Info = info;

			Grid = new Mesh();

			int pointCount = Info.HorizontalNodeCount * (Info.FirstLayerVerticalNodeCount + Info.SecondLayerVerticalNodeCount - 1);
			Points = new List<Point>();

			FB = new FirstBoundary();
		}

		public void Build()
		{
			BuildPoints();
			BuildTriangles();
			BuildBoundary();
		}

		void BuildPoints()
		{
			double HorizontalStep = Info.Width / (Info.HorizontalNodeCount - 1);
			double FirstLayerVerticalStep = Info.FirstLayerHeight / (Info.FirstLayerVerticalNodeCount - 1);
			double SecondLayerVerticalStep = Info.SecondLayerHeight / (Info.SecondLayerVerticalNodeCount - 1);

			for (int i = 0; i < Info.FirstLayerVerticalNodeCount; i++)
			{
				for (int j = 0; j < Info.HorizontalNodeCount; j++)
				{
					Points.Add(new Point(Info.R0 + j * HorizontalStep, -(Info.Z0 + i * FirstLayerVerticalStep)));
				}
			}

			double Z1 = Info.Z0 + Info.FirstLayerHeight;

			for (int i = 0; i < Info.SecondLayerVerticalNodeCount - 1; i++)
			{
				for (int j = 0; j < Info.HorizontalNodeCount; j++)
				{
					Points.Add(new Point(Info.R0 + j * HorizontalStep, -(Z1 + (i + 1) * SecondLayerVerticalStep)));
				}
			}
		}

		void BuildTriangles()
		{
			for (int i = 0; i < Info.FirstLayerVerticalNodeCount - 1; i++)
			{
				for (int j = 0; j < Info.HorizontalNodeCount - 1; j++)
				{
					FiniteElement e1 = new FiniteElement();
					FiniteElement e2 = new FiniteElement();

					int p1 = i * Info.HorizontalNodeCount + j;
					int p2 = i * Info.HorizontalNodeCount + j + 1;
					int p3 = (i + 1) * Info.HorizontalNodeCount + j;
					int p4 = (i + 1) * Info.HorizontalNodeCount + j + 1;

					e1.Vertices[0] = p1;
					e1.Vertices[1] = p2;
					e1.Vertices[2] = p3;
					e1.Material = Info.Materials[0];

					e2.Vertices[0] = p2;
					e2.Vertices[1] = p4;
					e2.Vertices[2] = p3;
					e2.Material = Info.Materials[0];

					Grid.Elements.Add(e1);
					Grid.Elements.Add(e2);
				}
			}

			for (int i = Info.FirstLayerVerticalNodeCount - 1; i < Info.FirstLayerVerticalNodeCount + Info.SecondLayerVerticalNodeCount - 2; i++)
			{
				for (int j = 0; j < Info.HorizontalNodeCount - 1; j++)
				{
					FiniteElement e1 = new FiniteElement();
					FiniteElement e2 = new FiniteElement();

					int p1 = i * Info.HorizontalNodeCount + j;
					int p2 = i * Info.HorizontalNodeCount + j + 1;
					int p3 = (i + 1) * Info.HorizontalNodeCount + j;
					int p4 = (i + 1) * Info.HorizontalNodeCount + j + 1;

					e1.Vertices[0] = p1;
					e1.Vertices[1] = p2;
					e1.Vertices[2] = p3;
					e1.Material = Info.Materials[1];

					e2.Vertices[0] = p2;
					e2.Vertices[1] = p4;
					e2.Vertices[2] = p3;
					e2.Material = Info.Materials[1];

					Grid.Elements.Add(e1);
					Grid.Elements.Add(e2);
				}
			}
		}

		void BuildBoundary()
		{
			// Top
			if (Info.Conditions[AreaSide.Top].Item1 == ConditionType.First)
			{
				for (int i = 0; i < Info.HorizontalNodeCount - 1; i++)
				{
					Edge edge = new Edge();
					edge.V1 = i;
					edge.V4 = i + 1;
					edge.Function = Info.Conditions[AreaSide.Top].Item2;
					FB.Edges.Add(edge);
				}
			}

			// Bottom
			if (Info.Conditions[AreaSide.Bottom].Item1 == ConditionType.First)
			{
				for (int i = 0; i < Info.HorizontalNodeCount - 1; i++)
				{
					Edge edge = new Edge();
					edge.V1 = Info.HorizontalNodeCount * (Info.FirstLayerVerticalNodeCount + Info.SecondLayerVerticalNodeCount - 2) + i;
					edge.V4 = Info.HorizontalNodeCount * (Info.FirstLayerVerticalNodeCount + Info.SecondLayerVerticalNodeCount - 2) + i + 1;
					edge.Function = Info.Conditions[AreaSide.Top].Item2;
					FB.Edges.Add(edge);
				}
			}

			// Left
			if (Info.Conditions[AreaSide.Left].Item1 == ConditionType.First)
			{
				for (int i = 0; i < Info.FirstLayerVerticalNodeCount + Info.SecondLayerVerticalNodeCount - 2; i++)
				{
					Edge edge = new Edge();
					edge.V1 = i * Info.HorizontalNodeCount;
					edge.V4 = (i + 1) * Info.HorizontalNodeCount;
					edge.Function = Info.Conditions[AreaSide.Top].Item2;
					FB.Edges.Add(edge);
				}
			}

			//Right
			if (Info.Conditions[AreaSide.Right].Item1 == ConditionType.First)
			{
				for (int i = 0; i < Info.FirstLayerVerticalNodeCount + Info.SecondLayerVerticalNodeCount - 2; i++)
				{
					Edge edge = new Edge();
					edge.V1 = (i + 1) * Info.HorizontalNodeCount - 1;
					edge.V4 = (i + 2) * Info.HorizontalNodeCount - 1;
					edge.Function = Info.Conditions[AreaSide.Top].Item2;
					FB.Edges.Add(edge);
				}
			}
		}
	}
}
