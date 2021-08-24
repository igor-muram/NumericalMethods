#pragma once

#include <vector>

using namespace std;

class Matrix
{
public:
	void Multiply(vector<double>& vec, vector<double>& res)
	{
		for (int i = 0; i < N; i++)
		{
			res[i] = vec[i] * DI[i];
			for (int k = IA[i]; k < IA[i + 1]; k++)
			{
				int j = JA[k];
				res[i] += AL[k] * vec[j];
				res[j] += AU[k] * vec[i];
			}
		}
	}

	void Forward(vector<double>& x, vector<double>& b)
	{
		for (int i = 0; i < N; i++)
		{
			double sum = 0;
			int i0 = IA[i], i1 = IA[i + 1];
			for (int k = i0; k < i1; k++)
			{
				int j = JA[k];
				sum += AL[k] * x[j];
			}
			x[i] = (b[i] - sum) / DI[i];
		}
	}

	void Backward(vector<double>& x, vector<double>& b)
	{
		for (int i = 0; i < N; i++)
			x[i] = b[i];

		for (int i = N - 1; i >= 0; i--)
		{
			int i0 = IA[i], i1 = IA[i + 1];
			for (int k = i0; k < i1; k++)
			{
				int j = JA[k];
				x[j] -= AU[k] * x[i];

			}
		}
	}

	void LUFactorization(Matrix& LU)
	{

		LU.N = N;
		LU.IA = IA;
		LU.JA = JA;

		LU.AL.resize(IA[N]);
		LU.AU.resize(IA[N]);
		LU.DI.resize(N);


		for (int i = 0; i < N; i++)
		{
			double sumD = 0;
			int i0 = IA[i], i1 = IA[i + 1];

			for (int k = i0; k < i1; k++)
			{
				double sumL = 0, sumU = 0;
				int j = JA[k];

				// Calculate L[i][j], U[j][i]
				int j0 = IA[j], j1 = IA[j + 1];

				int kl = i0, ku = j0;

				while(kl < i1 && ku < j1)
				{
					int j_kl = JA[kl];
					int j_ku = JA[ku];

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

				LU.AL[k] = AL[k] - sumL;
				LU.AU[k] = AU[k] - sumU;
				LU.AU[k] /= LU.DI[j];

				// Calculate sum for DI[i]
				sumD += LU.AL[k] * LU.AU[k];
			}

			// Calculate DI[i]
			LU.DI[i] = DI[i] - sumD;
		}

	}

public:
	int N;
	vector<double> DI, AL, AU;
	vector<int> IA, JA;
};

