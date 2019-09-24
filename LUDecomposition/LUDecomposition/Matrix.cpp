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

void Compute(int& N, real* DI, real* AL, real* AU, int* IA, real* b)
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