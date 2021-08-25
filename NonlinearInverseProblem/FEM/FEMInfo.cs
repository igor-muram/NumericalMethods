using System;
using System.Collections.Generic;

namespace FEM
{
	public class FEMInfo
	{
		public const int BasisSize = 10;

		public struct DerComp
		{
			public DerComp(int gradNo, double coefficient, int v1, int v2, int v3)
			{
				GradNo = gradNo;
				Coefficient = coefficient;
				V1 = v1;
				V2 = v2;
				V3 = v3;
			}
			public int GradNo { get; set; }
			public double Coefficient { get; set; }
			public int V1 { get; set; }
			public int V2 { get; set; }
			public int V3 { get; set; }
		};

		public struct PsiComp
		{
			public PsiComp(double coefficient, int v1, int v2, int v3)
			{
				Coefficient = coefficient;
				V1 = v1;
				V2 = v2;
				V3 = v3;
			}

			public double Coefficient { get; set; }
			public int V1 { get; set; }
			public int V2 { get; set; }
			public int V3 { get; set; }
		};

		public struct LocalCompG
		{
			public LocalCompG(int grad1, int grad2, double coefficient, int rnum)
			{
				Grad1 = grad1;
				Grad2 = grad2;
				Coefficient = coefficient;
				Rnum = rnum;
			}

			public int Grad1 { get; set; }
			public int Grad2 { get; set; }
			public double Coefficient { get; set; }
			public int Rnum { get; set; }
		};

		public struct LocalCompM
		{
			public LocalCompM(double coefficient, int rnum)
			{
				Coefficient = coefficient;
				Rnum = rnum;
			}

			public double Coefficient { get; set; }
			public int Rnum { get; set; }
		};

		static public List<List<DerComp>> Ders = new List<List<DerComp>>()
		{
			new List<DerComp>() { new DerComp(0, 13.5, 2, 0, 0),  new DerComp(0, -9.0, 1, 0, 0),   new DerComp(0, 1.0, 0, 0, 0) },
			new List<DerComp>() { new DerComp(1, 13.5, 0, 2, 0),  new DerComp(1, -9.0, 0, 1, 0),   new DerComp(1, 1.0, 0, 0, 0) },
			new List<DerComp>() { new DerComp(2, 13.5, 0, 0, 2),  new DerComp(2, -9.0, 0, 0, 1),   new DerComp(2, 1.0, 0, 0, 0) },

			new List<DerComp>() { new DerComp(0, 27.0, 1, 1, 0),  new DerComp(0, -4.5, 0, 1, 0),   new DerComp(1, 13.5, 2, 0, 0),   new DerComp(1, -4.5, 1, 0, 0) }, // 4
			new List<DerComp>() { new DerComp(0, 13.5, 0, 2, 0),  new DerComp(0, -4.5, 0, 1, 0),   new DerComp(1, 27.0, 1, 1, 0),   new DerComp(1, -4.5, 1, 0, 0) }, // 5

			new List<DerComp>() { new DerComp(0, 27.0, 1, 0, 1),  new DerComp(0, -4.5, 0, 0, 1),   new DerComp(2, 13.5, 2, 0, 0),   new DerComp(2, -4.5, 1, 0, 0) }, // 9
			new List<DerComp>() { new DerComp(0, 13.5, 0, 0, 2),  new DerComp(0, -4.5, 0, 0, 1),   new DerComp(2, 27.0, 1, 0, 1),   new DerComp(2, -4.5, 1, 0, 0) }, // 8

			new List<DerComp>() { new DerComp(1, 27.0, 0, 1, 1),  new DerComp(1, -4.5, 0, 0, 1),   new DerComp(2, 13.5, 0, 2, 0),   new DerComp(2, -4.5, 0, 1, 0) }, // 6
			new List<DerComp>() { new DerComp(1, 13.5, 0, 0, 2),  new DerComp(1, -4.5, 0, 0, 1),   new DerComp(2, 27.0, 0, 1, 1),   new DerComp(2, -4.5, 0, 1, 0) }, // 7
			new List<DerComp>() { new DerComp(0, 27.0, 0, 1, 1),  new DerComp(1, 27.0, 1, 0, 1),   new DerComp(2, 27.0, 1, 1, 0) }
		};

		static public List<List<PsiComp>> Basis = new List<List<PsiComp>>()
		{
			new List<PsiComp>() { new PsiComp(4.5, 3, 0, 0),      new PsiComp(-4.5, 2, 0, 0),   new PsiComp(1.0, 1, 0, 0) },
			new List<PsiComp>() { new PsiComp(4.5, 0, 3, 0),      new PsiComp(-4.5, 0, 2, 0),   new PsiComp(1.0, 0, 1, 0) },
			new List<PsiComp>() { new PsiComp(4.5, 0, 0, 3),      new PsiComp(-4.5, 0, 0, 2),   new PsiComp(1.0, 0, 0, 1) },

			new List<PsiComp>() { new PsiComp(13.5, 2, 1, 0),     new PsiComp(-4.5, 1, 1, 0) },
			new List<PsiComp>() { new PsiComp(13.5, 1, 2, 0),     new PsiComp(-4.5, 1, 1, 0) },

			new List<PsiComp>() { new PsiComp(13.5, 2, 0, 1),     new PsiComp(-4.5, 1, 0, 1) },
			new List<PsiComp>() { new PsiComp(13.5, 1, 0, 2),     new PsiComp(-4.5, 1, 0, 1) },

			new List<PsiComp>() { new PsiComp(13.5, 0, 2, 1),     new PsiComp(-4.5, 0, 1, 1) },
			new List<PsiComp>() { new PsiComp(13.5, 0, 1, 2),     new PsiComp(-4.5, 0, 1, 1) },
			new List<PsiComp>() { new PsiComp(27.0, 1, 1, 1) }
		};

		static public Func<double, double, double, double>[] LBasis = new Func<double, double, double, double>[BasisSize]
		{
			(double L1, double L2, double L3) => 0.5 * L1 * (3 * L1 - 1) * (3 * L1 - 2),
			(double L1, double L2, double L3) => 0.5 * L2 * (3 * L2 - 1) * (3 * L2 - 2),
			(double L1, double L2, double L3) => 0.5 * L3 * (3 * L3 - 1) * (3 * L3 - 2),
			(double L1, double L2, double L3) => 9.0 * L1 * L2 * (3 * L1 - 1) / 2.0,
			(double L1, double L2, double L3) => 9.0 * L1 * L2 * (3 * L2 - 1) / 2.0,
			(double L1, double L2, double L3) => 9.0 * L2 * L3 * (3 * L2 - 1) / 2.0,
			(double L1, double L2, double L3) => 9.0 * L2 * L3 * (3 * L3 - 1) / 2.0,
			(double L1, double L2, double L3) => 9.0 * L3 * L1 * (3 * L3 - 1) / 2.0,
			(double L1, double L2, double L3) => 9.0 * L3 * L1 * (3 * L1 - 1) / 2.0,
			(double L1, double L2, double L3) => 27.0 * L1 * L2 * L3
		};
	}
}
