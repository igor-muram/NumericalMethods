using ElectromagneticProblem.Matrix;
using MathUtility;
using System;

namespace ElectromagneticProblem.Solvers
{
   public enum SolverTypes { LOSLU, LOSLLT, CGMLLT }

   public interface ISolver
   {
      public SolverTypes SolverType { get; }
      public int MaxIterCount { get; set; }
      public int IterCount { get; set; }
      public double Eps { get; set; }
      public double Difference { get; set; }

      public double[] Solve(IMatrix matrix, double[] b);
   }

   public class LOSLLT : ISolver
   {
      public SolverTypes SolverType => SolverTypes.LOSLLT;

      public int MaxIterCount { get; set; } = 100000;
      public int IterCount { get; set; } = 0;

      public double Eps { get; set; } = 1.0e-15;
      public double Difference { get; set; } = 0.0;

      int N { get; set; } = 0;

      double[] Ax { get; set; } = null;
      double[] r { get; set; } = null;
      double[] z { get; set; } = null;
      double[] p { get; set; } = null;
      double[] temp { get; set; } = null;
      double[] xPrev { get; set; } = null;

      RawMatrix LLT { get; set; }
      struct RawMatrix
      {
         public int N { get; set; }
         public double[] DI { get; set; }
         public double[] AL { get; set; }
         public int[] IA { get; set; }
         public int[] JA { get; set; }
      }

      public double[] Solve(IMatrix matrix, double[] B)
      {
         N = matrix.N;
         InitAuxVectors(N);
         LLTFactorization(matrix);

         double[] x = new double[N];

         // Calculate r0
         matrix.Multiply(x, Ax);
         for (int i = 0; i < N; i++)
            r[i] = B[i] - Ax[i];

         Forward(LLT, r, r);

         //Calculate z0
         Backward(LLT, z, r);

         // Calculate p0
         matrix.Multiply(z, p);
         Forward(LLT, p, p);

         Difference = VectorUtils.DotProduct(r, r);

         while (IterCount < MaxIterCount && Difference >= Eps && VectorUtils.Error(x, xPrev) >= 1.0e-10)
         {
            // Calculate alpha
            double dotP = VectorUtils.DotProduct(p, p);
            double a = VectorUtils.DotProduct(p, r) / dotP;

            // Calculate xk, rk
            for (int i = 0; i < N; i++)
            {
               xPrev[i] = x[i];
               x[i] += a * z[i];
               r[i] -= a * p[i];
            }

            // Calculate beta
            Backward(LLT, Ax, r);
            matrix.Multiply(Ax, temp);
            Forward(LLT, Ax, temp);
            double b = -VectorUtils.DotProduct(p, Ax) / dotP;

            // Calculate zk, pk
            Backward(LLT, temp, r);
            for (int i = 0; i < N; i++)
            {
               z[i] = temp[i] + b * z[i];
               p[i] = Ax[i] + b * p[i];
            }

            // Calculate difference
            Difference = VectorUtils.DotProduct(r, r);

            IterCount++;
         }

         return x;
      }

      void LLTFactorization(IMatrix matrix)
      {
         RawMatrix LLT = new RawMatrix();
         LLT.N = matrix.N;
         LLT.IA = new int[LLT.N + 1];

         for (int i = 0; i < matrix.N + 1; i++)
            LLT.IA[i] = matrix.IA[i];

         LLT.AL = new double[LLT.IA[LLT.N]];
         LLT.JA = new int[LLT.IA[LLT.N]];
         LLT.DI = new double[LLT.N];

         for (int i = 0; i < matrix.IA[matrix.N]; i++)
            LLT.JA[i] = matrix.JA[i];

         for (int i = 0; i < matrix.N; i++)
         {
            double sumD = 0;
            int i0 = matrix.IA[i], i1 = matrix.IA[i + 1];

            for (int k = i0; k < i1; k++)
            {
               double sumL = 0, sumU = 0;
               int j = matrix.JA[k];

               // Calculate L[i][j], U[j][i]
               int j0 = matrix.IA[j], j1 = matrix.IA[j + 1];

               int kl = i0, ku = j0;

               for (; kl < i1 && ku < j1;)
               {
                  int j_kl = matrix.JA[kl];
                  int j_ku = matrix.JA[ku];

                  if (j_kl == j_ku)
                  {
                     sumL += LLT.AL[kl] * LLT.AL[ku];
                     kl++;
                     ku++;
                  }
                  if (j_kl > j_ku)
                     ku++;
                  if (j_kl < j_ku)
                     kl++;
               }

               LLT.AL[k] = (matrix.AL[k] - sumL) / LLT.DI[j];

               // Calculate sum for DI[i]
               sumD += LLT.AL[k] * LLT.AL[k];
            }

            // Calculate DI[i]
            LLT.DI[i] = Math.Sqrt(matrix.DI[i] - sumD);
         }

         this.LLT = LLT;
      }

      void InitAuxVectors(int N)
      {
         Ax = new double[N];
         r = new double[N];
         z = new double[N];
         p = new double[N];
         temp = new double[N];
         xPrev = new double[N];

         for (int i = 0; i < N; i++)
            xPrev[i] = 1.0;
      }

      void Forward(RawMatrix A, double[] x, double[] b)
      {
         double[] di = A.DI;
         double[] al = A.AL;
         int[] ia = A.IA;
         int[] ja = A.JA;
         int N = A.N;


         for (int i = 0; i < N; i++)
         {
            double sum = 0;
            int i0 = ia[i], i1 = ia[i + 1];
            for (int k = i0; k < i1; k++)
            {
               int j = ja[k];
               sum += al[k] * x[j];
            }
            x[i] = (b[i] - sum) / di[i];
         }
      }

      void Backward(RawMatrix A, double[] x, double[] b)
      {
         double[] di = A.DI;
         double[] al = A.AL;
         int[] ia = A.IA;
         int[] ja = A.JA;
         int N = A.N;

         for (int i = 0; i < N; i++)
            x[i] = b[i];

         for (int i = N - 1; i >= 0; i--)
         {
            int i0 = ia[i], i1 = ia[i + 1];
            x[i] /= di[i];
            for (int k = i0; k < i1; k++)
            {
               int j = ja[k];
               x[j] -= al[k] * x[i];

            }
         }
      }
   }

   public class LOSLU : ISolver
   {
      public SolverTypes SolverType => SolverTypes.LOSLU;

      public int MaxIterCount { get; set; } = 100000;
      public int IterCount { get; set; } = 0;

      public double Eps { get; set; } = 1.0e-15;
      public double Difference { get; set; } = 0.0;

      int N { get; set; } = 0;

      double[] Ax { get; set; } = null;
      double[] r { get; set; } = null;
      double[] z { get; set; } = null;
      double[] p { get; set; } = null;
      double[] temp { get; set; } = null;
      double[] xPrev { get; set; } = null;

      RawMatrix LU { get; set; }

      struct RawMatrix
      {
         public int N { get; set; }
         public double[] DI { get; set; }
         public double[] AL { get; set; }
         public double[] AU { get; set; }
         public int[] IA { get; set; }
         public int[] JA { get; set; }
      }

      public double[] Solve(IMatrix matrix, double[] B)
      {
         N = matrix.N;
         InitAuxVectors(N);
         LUFactorization(matrix);

         double[] x = new double[N];

         // Calculate r0
         matrix.Multiply(x, Ax);
         for (int i = 0; i < N; i++)
            r[i] = B[i] - Ax[i];

         Forward(LU, r, r);

         //Calculate z0
         Backward(LU, z, r);

         // Calculate p0
         matrix.Multiply(z, p);
         Forward(LU, p, p);

         Difference = VectorUtils.DotProduct(r, r);

         while (IterCount < MaxIterCount && Difference >= Eps && VectorUtils.Error(x, xPrev) >= 1.0e-10)
         {
            // Calculate alpha
            double dotP = VectorUtils.DotProduct(p, p);
            double a = VectorUtils.DotProduct(p, r) / dotP;

            // Calculate xk, rk
            for (int i = 0; i < N; i++)
            {
               xPrev[i] = x[i];
               x[i] += a * z[i];
               r[i] -= a * p[i];
            }

            // Calculate beta
            Backward(LU, Ax, r);
            matrix.Multiply(Ax, temp);
            Forward(LU, Ax, temp);
            double b = -VectorUtils.DotProduct(p, Ax) / dotP;

            // Calculate zk, pk
            Backward(LU, temp, r);
            for (int i = 0; i < N; i++)
            {
               z[i] = temp[i] + b * z[i];
               p[i] = Ax[i] + b * p[i];
            }

            // Calculate difference
            Difference = VectorUtils.DotProduct(r, r);

            IterCount++;
         }

         return x;
      }

      void LUFactorization(IMatrix matrix)
      {
         RawMatrix LU = new RawMatrix();
         LU.N = matrix.N;
         LU.IA = new int[LU.N + 1];

         for (int i = 0; i < matrix.N + 1; i++)
            LU.IA[i] = matrix.IA[i];

         LU.AL = new double[LU.IA[LU.N]];
         LU.AU = new double[LU.IA[LU.N]];
         LU.JA = new int[LU.IA[LU.N]];
         LU.DI = new double[LU.N];

         for (int i = 0; i < matrix.IA[matrix.N]; i++)
            LU.JA[i] = matrix.JA[i];

         for (int i = 0; i < matrix.N; i++)
         {
            double sumD = 0;
            int i0 = matrix.IA[i], i1 = matrix.IA[i + 1];

            for (int k = i0; k < i1; k++)
            {
               double sumL = 0, sumU = 0;
               int j = matrix.JA[k];

               // Calculate L[i][j], U[j][i]
               int j0 = matrix.IA[j], j1 = matrix.IA[j + 1];

               int kl = i0, ku = j0;

               for (; kl < i1 && ku < j1;)
               {
                  int j_kl = matrix.JA[kl];
                  int j_ku = matrix.JA[ku];

                  if (j_kl == j_ku)
                  {
                     sumL += LU.AL[kl] * LU.AU[ku];
                     sumU += LU.AU[kl] * LU.AL[ku];
                     kl++;
                     ku++;
                  }
                  if (j_kl > j_ku)
                     ku++;
                  if (j_kl < j_ku)
                     kl++;
               }

               LU.AL[k] = matrix.AL[k] - sumL;
               LU.AU[k] = matrix.AU[k] - sumU;
               LU.AU[k] /= LU.DI[j];

               // Calculate sum for DI[i]
               sumD += LU.AL[k] * LU.AU[k];
            }

            // Calculate DI[i]
            LU.DI[i] = matrix.DI[i] - sumD;
         }

         this.LU = LU;
      }

      void InitAuxVectors(int N)
      {
         Ax = new double[N];
         r = new double[N];
         z = new double[N];
         p = new double[N];
         temp = new double[N];
         xPrev = new double[N];

         for (int i = 0; i < N; i++)
            xPrev[i] = 1.0;
      }

      void Forward(RawMatrix A, double[] x, double[] b)
      {
         double[] di = A.DI;
         double[] al = A.AL;
         double[] au = A.AU;
         int[] ia = A.IA;
         int[] ja = A.JA;
         int N = A.N;


         for (int i = 0; i < N; i++)
         {
            double sum = 0;
            int i0 = ia[i], i1 = ia[i + 1];
            for (int k = i0; k < i1; k++)
            {
               int j = ja[k];
               sum += al[k] * x[j];
            }
            x[i] = (b[i] - sum) / di[i];
         }
      }

      void Backward(RawMatrix A, double[] x, double[] b)
      {
         double[] di = A.DI;
         double[] al = A.AL;
         double[] au = A.AU;
         int[] ia = A.IA;
         int[] ja = A.JA;
         int N = A.N;

         for (int i = 0; i < N; i++)
            x[i] = b[i];

         for (int i = N - 1; i >= 0; i--)
         {
            int i0 = ia[i], i1 = ia[i + 1];
            for (int k = i0; k < i1; k++)
            {
               int j = ja[k];
               x[j] -= au[k] * x[i];

            }
         }
      }
   }


    public class CGMLLT : ISolver
    {
        public SolverTypes SolverType => SolverTypes.CGMLLT;

        public int MaxIterCount { get; set; } = 100000;
        public int IterCount { get; set; } = 0;
        public double Eps { get; set; } = 1.0e-15;
        public double Difference { get; set; } = 0.0;

        int N = 0;

        double[] Ax = null;
        double[] r = null;
        double[] z = null;
        double[] xPrev = null;

        RawMatrix LLT { get; set; }
        struct RawMatrix
        {
            public int N { get; set; }
            public double[] DI { get; set; }
            public double[] AL { get; set; }
            public int[] IA { get; set; }
            public int[] JA { get; set; }
        }

        public double[] Solve(IMatrix matrix, double[] B)
        {
            N = matrix.N;
            InitAuxVectors(N);
            LLTFactorization(matrix);

            double[] x = new double[N];

            matrix.Multiply(x, Ax);
            for (int i = 0; i < N; i++)
                r[i] = B[i] - Ax[i];

            Forward(LLT, z, r);
            Backward(LLT, z, z);

            Difference = VectorUtils.Norm(r) / VectorUtils.Norm(B);

            double dot1 = 0;
            double dot2 = 0;

            dot1 = VectorUtils.DotProduct(z, r);

            while (IterCount < MaxIterCount && Difference >= Eps && VectorUtils.Error(x, xPrev) > 1.0e-10)
            {
                matrix.Multiply(z, Ax);

                double a = dot1 / VectorUtils.DotProduct(Ax, z);

                for (int i = 0; i < N; i++)
                {
                    xPrev[i] = x[i];
                    x[i] += a * z[i];
                    r[i] -= a * Ax[i];
                }

                Forward(LLT, Ax, r);
                Backward(LLT, Ax, Ax);

                dot2 = VectorUtils.DotProduct(Ax, r);
                double b = dot2 / dot1;
                dot1 = dot2;

                for (int i = 0; i < N; i++)
                    z[i] = Ax[i] + b * z[i];

                Difference = VectorUtils.Norm(r) / VectorUtils.Norm(B);
                IterCount++;
            }

            return x;
        }

        void LLTFactorization(IMatrix matrix)
        {
            RawMatrix LLT = new RawMatrix();
            LLT.N = matrix.N;
            LLT.IA = new int[LLT.N + 1];

            for (int i = 0; i < matrix.N + 1; i++)
                LLT.IA[i] = matrix.IA[i];

            LLT.AL = new double[LLT.IA[LLT.N]];
            LLT.JA = new int[LLT.IA[LLT.N]];
            LLT.DI = new double[LLT.N];

            for (int i = 0; i < matrix.IA[matrix.N]; i++)
                LLT.JA[i] = matrix.JA[i];

            for (int i = 0; i < matrix.N; i++)
            {
                double sumD = 0;
                int i0 = matrix.IA[i], i1 = matrix.IA[i + 1];

                for (int k = i0; k < i1; k++)
                {
                    double sumL = 0, sumU = 0;
                    int j = matrix.JA[k];

                    // Calculate L[i][j], U[j][i]
                    int j0 = matrix.IA[j], j1 = matrix.IA[j + 1];

                    int kl = i0, ku = j0;

                    for (; kl < i1 && ku < j1;)
                    {
                        int j_kl = matrix.JA[kl];
                        int j_ku = matrix.JA[ku];

                        if (j_kl == j_ku)
                        {
                            sumL += LLT.AL[kl] * LLT.AL[ku];
                            kl++;
                            ku++;
                        }
                        if (j_kl > j_ku)
                            ku++;
                        if (j_kl < j_ku)
                            kl++;
                    }

                    LLT.AL[k] = (matrix.AL[k] - sumL) / LLT.DI[j];

                    // Calculate sum for DI[i]
                    sumD += LLT.AL[k] * LLT.AL[k];
                }

                // Calculate DI[i]
                LLT.DI[i] = Math.Sqrt(matrix.DI[i] - sumD);
            }

            this.LLT = LLT;
        }

        void InitAuxVectors(int N)
        {
            Ax = new double[N];
            r = new double[N];
            z = new double[N];
            xPrev = new double[N];

            for (int i = 0; i < N; i++)
                xPrev[i] = 1.0;
        }

        void Forward(RawMatrix A, double[] x, double[] b)
        {
            double[] di = A.DI;
            double[] al = A.AL;
            int[] ia = A.IA;
            int[] ja = A.JA;
            int N = A.N;


            for (int i = 0; i < N; i++)
            {
                double sum = 0;
                int i0 = ia[i], i1 = ia[i + 1];
                for (int k = i0; k < i1; k++)
                {
                    int j = ja[k];
                    sum += al[k] * x[j];
                }
                x[i] = (b[i] - sum) / di[i];
            }
        }

        void Backward(RawMatrix A, double[] x, double[] b)
        {
            double[] di = A.DI;
            double[] al = A.AL;
            int[] ia = A.IA;
            int[] ja = A.JA;
            int N = A.N;

            for (int i = 0; i < N; i++)
                x[i] = b[i];

            for (int i = N - 1; i >= 0; i--)
            {
                int i0 = ia[i], i1 = ia[i + 1];
                x[i] /= di[i];
                for (int k = i0; k < i1; k++)
                {
                    int j = ja[k];
                    x[j] -= al[k] * x[i];

                }
            }
        }
    }
}
