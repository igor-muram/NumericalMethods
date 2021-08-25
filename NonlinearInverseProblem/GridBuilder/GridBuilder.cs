using MathUtilities;
using System.Collections.Generic;

namespace GridBuilder
{
	public struct AreaInfo
	{
		public double R { get; set; }
		public double Z { get; set; }
		public double Width { get; set; }
		public double Height { get; set; }

		public double FirstLayerHeight { get; set; }
		public double SecondLayerHeight { get; set; }
	}

	public class GridBuilder
	{
		public AreaInfo Info { get; set; }

		public GridBuilder(AreaInfo info)
		{
			Info = info;
		}

		public void Build()
		{

		}
	}

	public class Grid
	{
		public List<Triangle> Triangles { get; set; }
		public List<Point> Points { get; set; }
	}
}
