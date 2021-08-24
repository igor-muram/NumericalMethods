using ElectromagneticProblem.Enviroment;
using ElectromagneticProblem.FEM;
using System;
using ElectromagneticProblem.Solvers;
using ElectromagneticProblem.Splines;

namespace ElectromagneticProblem
{
   class Program
   {
      static void Main(string[] args)
      {
         Area area = Area.FromFile(@"C:\repos\NumericalMethods\ElectromagneticProblem\ElectromagneticProblem\input\grid.txt");

         if (area != null)
         {
            Mesh mesh = new Mesh();
            area.BuildMesh(mesh);

            ProblemInfo info = new ProblemInfo();
            info.Mesh = mesh;
            info.SolverType = SolverTypes.LOSLLT;

            CubicHermiteSpline spline = CubicHermiteSpline.FromFile(@"C:\repos\NumericalMethods\ElectromagneticProblem\ElectromagneticProblem\input\mu.txt");

            NonlinearProblem problem = new NonlinearProblem(info, spline);
            problem.Solve();
            //Console.WriteLine("Az: " + problem.GetValueA(new Point(-1.53E-02, 3.50E-03)));
            //Console.WriteLine("|B|: " + problem.GetValueB(new Point(-1.53E-02, 3.50E-03)));
         }
      
         Console.ReadKey();
      }
   }
}
