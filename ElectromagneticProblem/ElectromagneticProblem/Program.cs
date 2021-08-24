using ElectromagneticProblem.Enviroment;
using ElectromagneticProblem.FEM;
using System;
using ElectromagneticProblem.Solvers;
using ElectromagneticProblem.Splines;
using MathUtility;

namespace ElectromagneticProblem
{
   class Program
   {
      static void Main(string[] args)
      {
         Area linearArea = Area.FromFile(@"C:\repos\NumericalMethods\ElectromagneticProblem\ElectromagneticProblem\input\linGrid.txt");
         Area nonlinearArea = Area.FromFile(@"C:\repos\NumericalMethods\ElectromagneticProblem\ElectromagneticProblem\input\nonlinGrid.txt");

         if (linearArea != null && nonlinearArea != null)
         {
            MuSpline spline = MuSpline.FromFile(@"C:\repos\NumericalMethods\ElectromagneticProblem\ElectromagneticProblem\input\mu.txt");
            nonlinearArea.SetSplinesForMaterials(spline);

            Mesh linearMesh = new Mesh();
            Mesh nonlinearMesh = new Mesh();
            linearArea.BuildMesh(linearMesh);
            nonlinearArea.BuildMesh(nonlinearMesh);

            NonlinearProblemInfo info = new NonlinearProblemInfo
            {
               LinearMesh = linearMesh,
               NonlinearMesh = nonlinearMesh,
               SolverType = SolverTypes.LOSLLT,
               Eps = 0.01,
               Delta = 1.0e-20,
               MaxIters = 100,
               DoOptimization = true 
            };

            NonlinearProblem problem = new NonlinearProblem(info);
            problem.Solve();

            Console.WriteLine();
            Console.WriteLine("Az: " + problem.GetValueA(new Point(-1.53E-02, 3.50E-03)));
            Console.WriteLine("Az: " + problem.GetValueA(new Point(0.0044, 0.0015)));
            Console.WriteLine();
            Console.WriteLine();
            Console.WriteLine();

            //info.DoOptimization = false;
            //SimpleIteration problem1 = new SimpleIteration(info);
            //problem1.Solve();

            //Console.WriteLine();
            //Console.WriteLine("Az: " + problem.GetValueA(new Point(-1.53E-02, 3.50E-03)));
            //Console.WriteLine();
            //Console.WriteLine("Az: " + problem.GetValueA(new Point(0.0044, 0.0015)));
         }
      
         Console.ReadKey();
      }
   }
}
