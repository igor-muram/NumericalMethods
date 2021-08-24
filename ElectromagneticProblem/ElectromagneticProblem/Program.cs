using ElectromagneticProblem.Enviroment;
using ElectromagneticProblem.FEM;
using System;
using System.IO;
using MathUtility;
using ElectromagneticProblem.Solvers;

namespace ElectromagneticProblem
{
   class Program
   {
      static void Main(string[] args)
      {
         ProblemInfo info = new ProblemInfo();
         info.AreaFilename = @"C:\repos\NumericalMethods\ElectromagneticProblem\ElectromagneticProblem\input\grid1.txt";
         info.SolverType = SolverTypes.LOSLLT;

         Problem problem = new Problem(info);
         problem.Solve();
         Console.WriteLine("Az: " + problem.GetValueA(new Point(-1.53E-02, 3.50E-03)));
         Console.WriteLine("|B|: " + problem.GetValueB(new Point(-1.53E-02, 3.50E-03)));

         Console.ReadKey();

         //CubicHermiteSpline spline = CubicHermiteSpline.FromFile(@"C:\Users\Igor\source\repos\ElectromagneticProblem\ElectromagneticProblem\input\mu.txt");

         //double a = 0.01;
         //double step = 0.00005;
         //double b = 4.9;

         //using (StreamWriter writer = new StreamWriter(@"C:\Users\Igor\source\repos\ElectromagneticProblem\ElectromagneticProblem\input\output.csv"))
         //{
         //	while (a < b)
         //	{
         //		writer.WriteLine($"{a}; {spline.GetValue(a)}");
         //		a += step;
         //	}
         //}
      }
   }
}
