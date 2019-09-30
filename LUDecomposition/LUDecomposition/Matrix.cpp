#include "Matrix.h"

using namespace std;

void InputSize(int& N, int& ALSize)
{
	ifstream in("input/Size.txt");
	in >> N >> ALSize;
	in.close();
}

void Input(int& N, int& ALSize, real *DI, real *AL, real *AU, real *b, int *IA)
{
	ifstream in("input/DI.txt");
	for (int i = 0; i < N; i++)
		in >> DI[i];
	in.close();

	in.open("input/AL.txt");
	for (int i = 0; i < ALSize; i++)
		in >> AL[i];
	in.close();

	in.open("input/AU.txt");
	for (int i = 0; i < ALSize; i++)
		in >> AU[i];
	in.close();

	in.open("input/IA.txt");
	for (int i = 0; i < N + 1; i++)
	{
		in >> IA[i];
		IA[i]--;
	}
	in.close();

	in.open("input/b.txt");
	for (int i = 0; i < N; i++)
		in >> b[i];
	in.close();
}

void BuildLU(int& N, real *DI, real *AL, real *AU, int *IA)
{
	for (int i = 0; i < N; i++)
	{
		realScal sumD = 0;										// Sum for calculating elements of array D
		int i0 = IA[i], i1 = IA[i + 1];
		int j = i - (i1 - i0);

		for (int k = i0; k < i1; k++, j++)
		{
			realScal sumL = 0, sumU = 0;						// Sums for calculating elements of arrays AL and AU

			// Calculation elements L[i][j] and U[j][i]
			int j0 = IA[j], j1 = IA[j + 1];
			int size_i = k - i0, size_j = j1 - j0;
			int diff = size_i - size_j;
			int kl = i0, ku = j0;

			(diff < 0) ? ku -= diff : kl += diff;

			for (; kl < k; kl++, ku++)
			{
				sumL += AL[kl] * AU[ku];
				sumU += AU[kl] * AL[ku];
			}

			AL[k] -= sumL;
			AU[k] -= sumU;
			AU[k] /= DI[j];

			// Accumulation of sum for DI[i]
			sumD += AL[k] * AU[k];
		}

		// Calculation DI[i]
		DI[i] -= sumD;
		sumD = 0;
	}
}

void Compute(int& N, real *DI, real *AL, real *AU, int *IA, real *b)
{
	real *y = b, *x = b;

	// Solution of Ly = b by direct bypass
	for (int i = 0; i < N; i++)
	{
		realScal sumL = 0;
		int j = i - 1;

		for (int k = IA[i + 1] - 1; k >= IA[i]; k--, j--)
			sumL += AL[k] * y[j];
		
		y[i] = (b[i] - sumL) / DI[i];
	}

	// Solution of Ux = y by reverse bypass
	for (int i = N - 1; i >= 0; i--)
	{
		int j = i - 1;
		x[i] = y[i];
		for (int k = IA[i + 1] - 1; k >= IA[i]; k--, j--)
			y[j] -= AU[k] * x[i];
	}
}

void Output(int& N, int& ALSize, real *DI, real *AL, real *AU, real *x)
{
	cout << "DI: ";
	for (int i = 0; i < N; i++)
		cout << DI[i] << ' ';
	cout << endl << endl;

	cout << "AL: ";
	for (int i = 0; i < ALSize; i++)
		cout << AL[i] << ' ';
	cout << endl << endl;

	cout << "AU: ";
	for (int i = 0; i < ALSize; i++)
		cout << AU[i] << ' ';
	cout << endl << endl;

	cout << "x: " << endl;
	for (int i = 0; i < N; i++)
		cout << scientific << x[i] << endl;
	cout << endl << endl;

	cout << "x* - x: " << endl;
	for (int i = 0; i < N; i++)
		cout << scientific << i + 1 - x[i] << endl;
	cout << endl << endl;
}

void HilbertMatrix(int& N, real *DI, real *AL, real *AU, int *IA)
{
	for (int i = 0; i < N; i++)
		DI[i] = 1.0 / (2 * i + 1);

	IA[0] = 0;
	for (int i = 0; i < N; i++)
		IA[i + 1] = IA[i] + i;

	for (int i = 1; i < N; i++)
		for (int j = 0; j < i; j++)
			AU[IA[i] + j] = AL[IA[i] + j] = 1.0 / (i + j + 1);
}

void ToTightFormat(int& N, real *DI, real *AL, real *AU, int *IA, real **A)
{
	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			A[i][j] = 0;

	for (int i = 0; i < N; i++)
		A[i][i] = DI[i];

	for (int i = 1; i < N; i++)
	{
		int i0 = IA[i], i1 = IA[i + 1];
		int j = i - (i1 - i0);

		for (int k = i0; k < i1; k++, j++)
		{
			A[i][j] = AL[k];
			A[j][i] = AU[k];
		}
	}
}

void Gauss(int& N, real *x, real *b, real **A)
{
	for (int k = 0; k < N - 1; k++)
	{
		// Search for the leading element
		real max = abs(A[k][k]);
		int m = k;
		for (int i = k + 1; i < N; i++)
			if (abs(A[k][i]) > max)
			{
				max = abs(A[k][i]);
				m = i;
			}

		// Swap b[m] and b[k]
		real temp = b[m];
		b[m] = b[k];
		b[k] = temp;

		// Swap k and m columns
		for (int j = k; j < N; j++)
		{
			temp = A[k][j];
			A[k][j] = A[m][j];
			A[m][j] = temp;
		}

		// Zeroing k column
		for (int i = k + 1; i < N; i++)
		{
			real t = A[i][k] / A[k][k];
			b[i] -= t * b[k];
			for (int j = k + 1; j < N; j++)
				A[i][j] -= t * A[k][j];
		}
	}

	// Calculation of x
	x[N - 1] = b[N - 1] / A[N - 1][N - 1];
	for (int k = N - 2; k >= 0; k--)
	{
		realScal sum = 0;
		for (int j = k + 1; j < N; j++)
			sum += A[k][j] * x[j];

		x[k] = (b[k] - sum) / A[k][k];
	}
}

void Multiply(int& N, real *DI, real *AL, real *AU, int *IA, real *vector, real *res)
{
	for (int i = 0; i < N; i++)
		res[i] = vector[i] * DI[i];

	for (int i = 1; i < N; i++)
	{
		int j = i - (IA[i + 1] - IA[i]);

		for (int k = IA[i]; k < IA[i + 1]; k++, j++)
		{
			res[i] += vector[j] * AL[k];
			res[j] += vector[i] * AU[k];
		}
	}
}