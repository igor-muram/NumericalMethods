#include "Solvers.h"
#include <cmath>

namespace Solvers
{
	struct AuxVectors
	{
		double* Ax = nullptr;
		double* r = nullptr;
		double* z = nullptr;
		double* p = nullptr;
		double* LU = nullptr;
		double* temp = nullptr;
		double* r0 = nullptr;
		double* y = nullptr;
		double* Ap = nullptr;
		double* Az = nullptr;
	};

	struct Matrix
	{
		int N;
		double* DI = nullptr;
		double* AL = nullptr;
		double* AU = nullptr;
		int* IA = nullptr;
		int* JA = nullptr;
	};

	int LOS_LU(Matrix& A, double* x, double* f, Matrix& LU, AuxVectors& aux, int maxiter, double eps, double& lastdiff);
	int BCG_LU(Matrix& A, double* x, double* f, Matrix& LU, AuxVectors& aux, int maxiter, double eps);
	void LUFactorization(Matrix& A, Matrix& LU);

	int LOS(SparseMatrix& A, std::vector<double>& x, std::vector<double>& b)
	{
		Matrix Raw;
		Raw.N = A.N;
		Raw.DI = A.DI.data();
		Raw.AL = A.AL.data();
		Raw.AU = A.AU.data();
		Raw.IA = A.IA.data();
		Raw.JA = A.JA.data();

		AuxVectors aux;

		aux.Ax = new double[A.N];
		aux.r = new double[A.N];
		aux.z = new double[A.N];
		aux.p = new double[A.N];
		aux.temp = new double[A.N];

		Matrix LU;
		LUFactorization(Raw, LU);
		std::vector<double> Braw = b;
		double lastdiff = 0;

		x.resize(A.N);
		return LOS_LU(Raw, x.data(), Braw.data(), LU, aux, 20000, 1.0e-7, lastdiff);

	}

	int BCG(SparseMatrix& A, std::vector<double>& x, std::vector<double>& b)
	{
		Matrix Raw;
		Raw.N = A.N;
		Raw.DI = A.DI.data();
		Raw.AL = A.AL.data();
		Raw.AU = A.AU.data();
		Raw.IA = A.IA.data();
		Raw.JA = A.JA.data();

		AuxVectors aux;

		aux.Ax = new double[A.N]{ 0.0 };
		aux.Ap = new double[A.N]{ 0.0 };
		aux.Az = new double[A.N]{ 0.0 };
		aux.r = new double[A.N]{ 0.0 };
		aux.z = new double[A.N]{ 0.0 };
		aux.p = new double[A.N]{ 0.0 };
		aux.y = new double[A.N]{ 0.0 };
		aux.r0 = new double[A.N]{ 0.0 };

		Matrix LU;
		LUFactorization(Raw, LU);
		double lastdiff = 0;

		x.resize(A.N);
		return BCG_LU(Raw, x.data(), b.data(), LU, aux, 20000, 1.0e-7);
	}

	void Multiply(Matrix& A, double* vec, double* res)
	{
		double* di = A.DI;
		double* al = A.AL;
		double* au = A.AU;
		int* ia = A.IA;
		int* ja = A.JA;
		int N = A.N;

		for (int i = 0; i < N; i++)
		{
			res[i] = vec[i] * di[i];
			for (int k = ia[i]; k < ia[i + 1]; k++)
			{
				int j = ja[k];
				res[i] += al[k] * vec[j];
				res[j] += au[k] * vec[i];
			}
		}
	}

	void MultiplyT(Matrix& A, double* vec, double* res)
	{
		double* di = A.DI;
		double* al = A.AU;
		double* au = A.AL;
		int* ia = A.IA;
		int* ja = A.JA;
		int N = A.N;

		for (int i = 0; i < N; i++)
		{
			res[i] = vec[i] * di[i];
			for (int k = ia[i]; k < ia[i + 1]; k++)
			{
				int j = ja[k];
				res[i] += al[k] * vec[j];
				res[j] += au[k] * vec[i];
			}
		}
	}

	void MultiplyU(Matrix& A, double* vec, double* res)
	{
		double* au = A.AU;
		int* ia = A.IA;
		int* ja = A.JA;
		int N = A.N;

		for (int i = 0; i < N; i++)
		{
			res[i] = vec[i];
			for (int k = ia[i]; k < ia[i + 1]; k++)
			{
				int j = ja[k];
				res[j] += au[k] * vec[i];
			}
		}
	}

	void LUFactorization(Matrix& A, Matrix& LU)
	{
		LU.N = A.N;
		LU.IA = new int[LU.N + 1];

		for (int i = 0; i < A.N + 1; i++)
			LU.IA[i] = A.IA[i];

		LU.AL = new double[LU.IA[LU.N]];
		LU.AU = new double[LU.IA[LU.N]];
		LU.JA = new int[LU.IA[LU.N]];
		LU.DI = new double[LU.N];

		for (int i = 0; i < A.IA[A.N]; i++)
			LU.JA[i] = A.JA[i];

		for (int i = 0; i < A.N; i++)
		{
			double sumD = 0;
			int i0 = A.IA[i], i1 = A.IA[i + 1];

			for (int k = i0; k < i1; k++)
			{
				double sumL = 0, sumU = 0;
				int j = A.JA[k];

				// Calculate L[i][j], U[j][i]
				int j0 = A.IA[j], j1 = A.IA[j + 1];

				int kl = i0, ku = j0;

				for (; kl < i1 && ku < j1; )
				{
					int j_kl = A.JA[kl];
					int j_ku = A.JA[ku];

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

				LU.AL[k] = A.AL[k] - sumL;
				LU.AU[k] = A.AU[k] - sumU;
				LU.AU[k] /= LU.DI[j];

				// Calculate sum for DI[i]
				sumD += LU.AL[k] * LU.AU[k];
			}

			// Calculate DI[i]
			LU.DI[i] = A.DI[i] - sumD;
		}

	}

	void Forward(Matrix& A, double* x, double* b, bool transposed)
	{
		double* di = A.DI;
		double* al = A.AL;
		double* au = A.AU;
		int* ia = A.IA;
		int* ja = A.JA;
		int N = A.N;

		if (transposed)
		{
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
		else
		{
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
	}

	void Backward(Matrix& A, double* x, double* b, bool transposed)
	{
		double* di = A.DI;
		double* al = A.AL;
		double* au = A.AU;
		int* ia = A.IA;
		int* ja = A.JA;
		int N = A.N;

		if (transposed)
		{
			for (int i = 0; i < N; i++)
			{
				double sum = 0;
				int i0 = ia[i], i1 = ia[i + 1];
				for (int k = i0; k < i1; k++)
				{
					int j = ja[k];
					sum += au[k] * x[j];
				}
				x[i] = b[i] - sum;
			}
		}
		else
		{
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

	double DotProduct(int N, double* a, double* b)
	{
		double sum = 0;
		for (int i = 0; i < N; i++)
			sum += a[i] * b[i];

		return sum;
	}

	int LOS_LU(Matrix& A, double* x, double* f, Matrix& LU, AuxVectors& aux, int maxiter, double eps, double& lastdiff)
	{
		double* Ax = aux.Ax;
		double* r = aux.r;
		double* z = aux.z;
		double* p = aux.p;
		double* temp = aux.temp;
		int N = A.N;

		// Calculate r0
		Multiply(A, x, Ax);
		for (int i = 0; i < N; i++)
			r[i] = f[i] - Ax[i];
		Forward(LU, r, r, false);

		//Calculate z0
		Backward(LU, z, r, false);

		// Calculate p0
		Multiply(A, z, p);
		Forward(LU, p, p, false);

		double diff = DotProduct(N, r, r);

		int k = 0;
		for (; k < maxiter && diff >= eps; k++)
		{
			// Calculate alpha
			double dotP = DotProduct(N, p, p);
			double a = DotProduct(N, p, r) / dotP;

			// Calculate xk, rk
			for (int i = 0; i < N; i++)
			{
				x[i] += a * z[i];
				r[i] -= a * p[i];
			}

			// Calculate beta
			Backward(LU, Ax, r, false);
			Multiply(A, Ax, temp);
			Forward(LU, Ax, temp, false);
			double b = -DotProduct(N, p, Ax) / dotP;

			// Calculate zk, pk
			Backward(LU, temp, r, false);
			for (int i = 0; i < N; i++)
			{
				z[i] = temp[i] + b * z[i];
				p[i] = Ax[i] + b * p[i];
			}

			// Calculate difference
			diff = DotProduct(N, r, r);
		}

		lastdiff = diff;
		return k;
	}

	int BCG_LU(Matrix& A, double* x, double* f, Matrix& LU, AuxVectors& aux, int maxiter, double eps)
	{
		double* r0 = aux.r0;
		double* r = aux.r;
		double* z = aux.z;
		double* p = aux.p;
		double* y = aux.y;
		double* Ax = aux.Ax;
		double* Ap = aux.Ap;
		double* Az = aux.Az;

		int N = A.N;

		// r(0)
		Multiply(A, x, Ax);
		for (int i = 0; i < N; i++)
			r0[i] = f[i] - Ax[i];
		Forward(LU, r0, r0, false);

		// z(0)
		Backward(LU, z, r0, false);

		// r = r(0)
		for (int i = 0; i < N; i++)
			r[i] = r0[i];

		// diff
		double diff = DotProduct(N, r0, r0);

		int k = 0;
		for (; k < maxiter && diff >= eps; k++)
		{
			// L(-1) * A * U(-1) * z(k - 1)
			Backward(LU, Ax, z, false);
			Multiply(A, Ax, Az);
			Forward(LU, Az, Az, false);

			// alpha(k)
			double dotP = DotProduct(N, r, r0);
			double a = dotP / DotProduct(N, Az, r0);

			// p(k)
			for (int i = 0; i < N; i++)
				p[i] = r[i] - a * Az[i];

			// L(-1) * A * U(-1) * p(k)
			Backward(LU, Ax, p, false);
			Multiply(A, Ax, Ap);
			Forward(LU, Ap, Ap, false);

			// gamma(k)
			double g = DotProduct(N, p, Ap) / DotProduct(N, Ap, Ap);

			// y(k)
			for (int i = 0; i < N; i++)
				y[i] = y[i] + a * z[i] + g * p[i];

			// r(k)
			for (int i = 0; i < N; i++)
				r[i] = p[i] - g * Ap[i];

			// beta(k)
			double b = a * DotProduct(N, r, r0) / (g * dotP);

			// z(k)
			for (int i = 0; i < N; i++)
				z[i] = r[i] + b * z[i] - b * g * Az[i];

			// diff
			diff = DotProduct(N, r, r);
		}

		// solution
		Backward(LU, x, y, false);
		return k;
	}
}