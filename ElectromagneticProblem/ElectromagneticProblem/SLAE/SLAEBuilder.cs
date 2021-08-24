using System;
using ElectromagneticProblem.Enviroment;
using ElectromagneticProblem.FEM;
using ElectromagneticProblem.Matrix;
using MathUtility;

namespace ElectromagneticProblem.SLAE
{
   class SLAEBuilder
   {
      public Mesh Mesh { get; set; } = null;

      double[,] localG1 = new double[4, 4]
         {
            {  2, -2,  1, -1 },
            { -2,  2, -1,  1 },
            {  1, -1,  2, -2 },
            { -1,  1, -2,  2 }
         };

      double[,] localG2 = new double[4, 4]
         {
            {  2,  1, -2, -1 },
            {  1,  2, -1, -2 },
            { -2, -1,  2,  1 },
            { -1, -2,  1,  2 }
         };

      double[,] localM = new double[4, 4]
         {
            { 4, 2, 2, 1 },
            { 2, 4, 1, 2 },
            { 2, 1, 4, 2 },
            { 1, 2, 2, 4 }
         };

      public SLAEBuilder(Mesh mesh)
      {
         Mesh = mesh;
      }

      public void Build(IMatrix A, double[] b)
      {
         foreach (FiniteElement e in Mesh.Elements)
         {
            double[,] localMatrix = BuildLocalMatrix(e);
            double[] localB = BuildLocalB(e);
            AddToGlobalMatrix(A, e, localMatrix);
            AddToGlobalB(b, e, localB);
         }

         AddFirstBoundary(A, b);
      }

      void AddFirstBoundary(IMatrix A, double[] b)
      {
         foreach (int index in Mesh.FirstBoundary.Indices)
         {
            A.Set(index, index, 1.0e+50);
            b[index] = Mesh.FirstBoundary.Value * 1.0e+50;
         }
      }

      double[,] BuildLocalMatrix(FiniteElement e)
      {
         Point p1 = Mesh.Points[e[0]];
         Point p2 = Mesh.Points[e[1]];
         Point p3 = Mesh.Points[e[2]];

         double hx = p2.X - p1.X;
         double hy = p3.Y - p1.Y;

         double[,] G = new double[FEMParameters.BasisSize, FEMParameters.BasisSize];

         double lambda = 1.0 / e.Material.mu;

         for (int i = 0; i < FEMParameters.BasisSize; i++)
            for (int j = 0; j < FEMParameters.BasisSize; j++)
               G[i, j] = lambda * (hy * localG1[i, j] / hx + hx * localG2[i, j] / hy) / 6;

         return G;
      }

      double[] BuildLocalB(FiniteElement e)
      {
         Point p1 = Mesh.Points[e[0]];
         Point p2 = Mesh.Points[e[1]];
         Point p3 = Mesh.Points[e[2]];

         double hx = p2.X - p1.X;
         double hy = p3.Y - p1.Y;

         double[] B = new double[FEMParameters.BasisSize];
         double J = e.Material.J;

         for (int i = 0; i < FEMParameters.BasisSize; i++)
            for (int j = 0; j < FEMParameters.BasisSize; j++)
               B[i] += localM[i, j] * J;

         for (int i = 0; i < FEMParameters.BasisSize; i++)
            B[i] *= hx * hy / 36;

         return B;
      }

      void AddToGlobalMatrix(IMatrix A, FiniteElement e, double[,] local)
      {
         for (int i = 0; i < FEMParameters.BasisSize; i++)
            for (int j = 0; j <= i; j++)
            {
               int iGlobal = e[i];
               int jGlobal = e[j];

               if (jGlobal > iGlobal)
                  (jGlobal, iGlobal) = (iGlobal, jGlobal);

               A.Add(iGlobal, jGlobal, local[i, j]);
            }
      }

      void AddToGlobalB(double[] B, FiniteElement e, double[] local)
      {
         for (int i = 0; i < FEMParameters.BasisSize; i++)
            B[e[i]] += local[i];
      }
   };
}
