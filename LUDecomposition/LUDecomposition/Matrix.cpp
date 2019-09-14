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
		real sumL = 0, sumU = 0, sumD = 0;					// Sums for calculating elements of arrays AL, AU and D
		int j = i - (IA[i + 1] - IA[i]);

		for (int k = IA[i]; k < IA[i + 1]; k++)
		{
			// Calculation L[i][j] and U[j][i]
			for (int m = 0; m < k - IA[i]; m++)
			{
				sumL += AL[IA[i] + m] * AU[IA[j] + m];
				sumU += AU[IA[i] + m] * AL[IA[j] + m];
			}

			AL[k] -= sumL;
			AU[k] -= sumU;
			AU[k] /= DI[j];

			// Accumulation of sum for DI[i]
			sumD += AL[k] * AU[k];

			sumL = sumU = 0;
			j++;
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
		real sumL = 0;
		int m = i - 1;

		for (int k = IA[i + 1] - 1; k >= IA[i]; k--)
		{
			sumL += AL[k] * y[m];
			m--;
		}
		
		y[i] = (b[i] - sumL) / DI[i];
	}

	// Solution of Ux = y by reverse bypass
	for (int i = N - 1; i > 0; i--)
	{
		int m = i - 1;
		for (int k = IA[i + 1] - 1; k >= IA[i]; k--)
		{
			x[m] = y[m] - AU[k] * x[i];
			m--;
		}
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

	cout << "x: ";
	for (int i = 0; i < N; i++)
		cout << x[i] << ' ';
	cout << endl << endl;
}