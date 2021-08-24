using System;
using System.Collections.Generic;
using System.IO;
using System.Globalization;
using System.Linq;
using MathUtility;
using ElectromagneticProblem.FEM;
using System.Collections;
using ElectromagneticProblem.Splines;

namespace ElectromagneticProblem.Enviroment
{
   public interface IMaterial
   {
      public Func<double, double> Mu { get; set; }
      public double J { get; set; }
   }

	public class Material : IMaterial
	{
      public Material(double mu, double J)
		{
         Mu = (double B) => mu;
         this.J = J;
		}

		public Func<double, double> Mu { get; set; }
		public double J { get; set; }
	}

	public class SplineMaterial : IMaterial
   {
      public SplineMaterial(ISpline spline, double J)
      {
         Mu = (double B) => spline.GetValue(B);
         this.J = J;
      }

		public Func<double, double> Mu { get; set; }
		public double J { get; set; }
	}

   public class Subarea
   {
      public double x1, x2;
      public double y1, y2;

      public IMaterial material;

      public static Subarea Parse(string data)
      {
         Subarea subarea = new Subarea();

         string[] tokens = data.Split(' ', StringSplitOptions.RemoveEmptyEntries);
         subarea.x1 = double.Parse(tokens[0], NumberStyles.Any, CultureInfo.InvariantCulture);
         subarea.x2 = double.Parse(tokens[1], NumberStyles.Any, CultureInfo.InvariantCulture);
         subarea.y1 = double.Parse(tokens[2], NumberStyles.Any, CultureInfo.InvariantCulture);
         subarea.y2 = double.Parse(tokens[3], NumberStyles.Any, CultureInfo.InvariantCulture);
         double mu = double.Parse(tokens[4], NumberStyles.Any, CultureInfo.InvariantCulture);
         double J = double.Parse(tokens[5], NumberStyles.Any, CultureInfo.InvariantCulture);
         int id = int.Parse(tokens[6], NumberStyles.Any, CultureInfo.InvariantCulture);

         mu *= 4 * Math.PI * 1.0e-7;

         subarea.material = new Material(mu, J);
         return subarea;
      }

      public bool Contains(Subarea subarea) => x1 <= subarea.x1 && subarea.x2 <= x2 && y1 <= subarea.y1 && subarea.y2 <= y2;
      public bool Contains(Interval x, Interval y) => x1 <= x.a && x.b <= x2 && y1 <= y.a && y.b <= y2;
   }

   public class Interval : IEnumerable
   {
      public double a;
      public double b;
      public double initialStep;
      public double q;
      public int dir;
      public int factor;
      public bool isUniform;

      public double[] Points { get; set; }

      public Interval(double a, double b, double initialStep, double q, int dir, int factor)
      {
         this.a = a;
         this.b = b;
         this.factor = factor;
         this.initialStep = initialStep;
         this.q = q;
         this.dir = dir;

         isUniform = q == 1.0;

         if (factor != 0)
         {
            this.initialStep /= Math.Pow(2.0, factor);

            if (!isUniform)
               this.q = Math.Pow(q, 1.0 / Math.Pow(2.0, factor));
			}


         Split();
      }

      void Split()
      {
         if (!isUniform)
         {
            int n = (int)Math.Ceiling(Math.Log((b - a) * (q - 1) / initialStep + 1) / Math.Log(q));
            initialStep = (b - a) * (q - 1) / (Math.Pow(q, n) - 1);

            List<double> steps = new List<double>();

            for (int i = 0; i < n; i++)
               steps.Add(initialStep * Math.Pow(q, i));

            Points = new double[steps.Count + 1];

            if (dir == 1)
            {
               double startPoint = a;
               Points[0] = startPoint;

               for (int i = 0; i < steps.Count; i++)
               {
                  Points[i + 1] = startPoint + steps[i];
                  startPoint += steps[i];
               }
            }
            else
            {
               double startPoint = a;
               Points[0] = startPoint;

               for (int i = steps.Count - 1; i >= 0; i--)
               {
                  Points[steps.Count - i] = startPoint + steps[i];
                  startPoint += steps[i];
               }
            }
         }
         else
         {
            int n = (int)Math.Round((b - a) / initialStep);
            initialStep = (b - a) / n;

            Points = new double[n + 1];
            for (int i = 0; i < n + 1; i++)
               Points[i] = a + i * initialStep;
         }
      }

      public IEnumerator GetEnumerator() => Points.GetEnumerator();
   }

   public class Area
   {
      public List<Subarea> Subareas { get; set; } = new List<Subarea>();
      Stack<List<Subarea>> Layers { get; set; } = new Stack<List<Subarea>>();
      public List<Interval> XIntervals { get; set; } = new List<Interval>();
      public List<Interval> YIntervals { get; set; } = new List<Interval>();
      public int XFactor { get; set; } = 0;
      public int YFactor { get; set; } = 0;

      public void BuildMesh(Mesh mesh)
      {
         List<double> xPoints = new List<double>();
         List<double> yPoints = new List<double>();

         foreach (var i in XIntervals)
            for (int k = 0; k < i.Points.Length - 1; k++)
               xPoints.Add(i.Points[k]);

         xPoints.Add(XIntervals.Last().Points.Last());

         foreach (var i in YIntervals)
            for (int k = 0; k < i.Points.Length - 1; k++)
               yPoints.Add(i.Points[k]);

         yPoints.Add(YIntervals.Last().Points.Last());

         mesh.NodeCount = xPoints.Count * yPoints.Count;

         // Set mesh points
         mesh.Points = new Point[xPoints.Count * yPoints.Count];
         int m = 0;
         foreach (double y in yPoints)
            foreach (double x in xPoints)
               mesh.Points[m++] = new Point(x, y);

         // Set mesh elements
         for (int j = 0; j < yPoints.Count - 1; j++)
         {
            for (int i = 0; i < xPoints.Count - 1; i++)
            {
               FiniteElement element = new FiniteElement(FEMParameters.ElementVerticesCount);

               element[0] = i + j * xPoints.Count;
               element[1] = i + 1 + j * xPoints.Count;
               element[2] = i + (j + 1) * xPoints.Count;
               element[3] = i + 1 + (j + 1) * xPoints.Count;

               mesh.Elements.Add(element);
            }
         }

         // Set materials for mesh elements
         foreach (var e in mesh.Elements)
         {
            double x1 = mesh.Points[e[0]].X;
            double x2 = mesh.Points[e[1]].X;
            double y1 = mesh.Points[e[0]].Y;
            double y2 = mesh.Points[e[2]].Y;

            Interval xInterval = null;
            foreach (var xInt in XIntervals)
               if (x1 >= xInt.a && x2 <= xInt.b)
               {
                  xInterval = xInt;
                  break;
               }

            Interval yInterval = null;
            foreach (var yInt in YIntervals)
               if (y1 >= yInt.a && y2 <= yInt.b)
               {
                  yInterval = yInt;
                  break;
               }

            if (xInterval == null || yInterval == null)
               throw new Exception("ERROR! Failed to find intervals for finite element");

            e.Material = FindSubarea(xInterval, yInterval).material;
         }

         mesh.FirstBoundary = new FirstNullBoundary();

         // Set mesh boundary conditions
         for (int i = 0; i < yPoints.Count; i++)
            mesh.FirstBoundary.Indices.Add(i * xPoints.Count);

         for (int i = 0; i < yPoints.Count; i++)
            mesh.FirstBoundary.Indices.Add(i * xPoints.Count + xPoints.Count - 1);

         for (int i = 1; i < xPoints.Count - 1; i++)
            mesh.FirstBoundary.Indices.Add(i + xPoints.Count * (yPoints.Count - 1));

         //for (int i = 1; i < xPoints.Count - 1; i++)
         //   mesh.FirstBoundary.Indices.Add(i);
      }

      void BuildLayers(Stack<List<Subarea>> layers, List<Subarea> subareas)
      {
         layers.Push(new List<Subarea>());

         foreach (var sub1 in subareas)
         {
            bool isContainted = false;
            foreach (var sub2 in subareas)
               if (sub1 != sub2)
                  if (sub2.Contains(sub1))
                  {
                     isContainted = true;
                     break;
                  }

            if (!isContainted)
               layers.Peek().Add(sub1);
         }

         foreach (var subarea in layers.Peek())
            subareas.Remove(subarea);

         if (subareas.Count != 0)
            BuildLayers(layers, subareas);
      }

      Subarea FindSubarea(Interval x, Interval y)
      {

         foreach (var layer in Layers)
         {
            foreach (var subarea in layer)
               if (subarea.Contains(x, y))
                  return subarea;
         }

         throw new Exception("ERROR! Couldn't find subarea for rectangle");
      }

      public static Area FromFile(string filename)
      {
         string[] data;
         try
         {
            data = File.ReadAllLines(filename);
         }
         catch
         {
            Console.WriteLine($"ERROR! Failed to open file {filename}");
            return null;
         }

         Area area = new Area();

         int rectCount = int.Parse(data[0], NumberStyles.Any, CultureInfo.InvariantCulture);

         for (int i = 0; i < rectCount; i++)
            area.Subareas.Add(Subarea.Parse(data[i + 1]));

         area.BuildLayers(area.Layers, area.Subareas.ToList());

         int intervalLines = 5;

         string[] factors = data.Last().Split(' ', StringSplitOptions.RemoveEmptyEntries);
         area.XFactor = int.Parse(factors[0], NumberStyles.Any);
         area.YFactor = int.Parse(factors[1], NumberStyles.Any);

         area.XIntervals = ParseIntervals(data.AsSpan(1 + rectCount + 1, intervalLines).ToArray(), area.XFactor);
         area.YIntervals = ParseIntervals(data.AsSpan(1 + rectCount + 1 + intervalLines + 1, intervalLines).ToArray(), area.YFactor);

         return area;
      }

      static List<Interval> ParseIntervals(string[] data, int factor)
      {
         List<Interval> intervals = new List<Interval>();
     
         string[] tokens = data[0].Split(' ', StringSplitOptions.RemoveEmptyEntries);

         double startPoint = double.Parse(tokens[0], NumberStyles.Any, CultureInfo.InvariantCulture);
         int intervalNum = int.Parse(tokens[1], NumberStyles.Any, CultureInfo.InvariantCulture);

         string[] points = data[1].Split(' ', StringSplitOptions.RemoveEmptyEntries);
         string[] initialSteps = data[2].Split(' ', StringSplitOptions.RemoveEmptyEntries);
         string[] Qs = data[3].Split(' ', StringSplitOptions.RemoveEmptyEntries);
         string[] dirs = data[4].Split(' ', StringSplitOptions.RemoveEmptyEntries);

         for (int i = 0; i < intervalNum; i++)
         {
            double endPoint = double.Parse(points[i], NumberStyles.Any, CultureInfo.InvariantCulture);
            double step = double.Parse(initialSteps[i], NumberStyles.Any, CultureInfo.InvariantCulture);
            double q = double.Parse(Qs[i], NumberStyles.Any, CultureInfo.InvariantCulture);
            int dir = int.Parse(dirs[i], NumberStyles.Any, CultureInfo.InvariantCulture);

            intervals.Add(new Interval(startPoint, endPoint, step, q, dir, factor));
            startPoint = endPoint;
         }

         return intervals;
      }

      Area() { }
   }
}
