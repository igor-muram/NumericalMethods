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
		real sumL = 0, sumU = 0, sumD = 0;					// Суммы для вычисления элементов массивов AL, AU и D
		int j = i - (IA[i + 1] - IA[i]);

		for (int k = IA[i]; k < IA[i + 1]; k++)
		{
			// Вычисление L[i][j] и U[j][i]
			for (int m = 0; m < k - IA[i]; m++)
			{
				sumL += AL[IA[i] + m] * AU[IA[j] + m];
				sumU += AU[IA[i] + m] * AL[IA[j] + m];
			}

			AL[k] -= sumL;
			AU[k] -= sumU;
			AU[k] /= DI[j];

			// Накопление суммы для DI[i]
			sumD += AL[k] * AU[k];

			sumL = sumU = 0;
			j++;
		}

		// Вычисление DI[i]
		DI[i] -= sumD;
		sumD = 0;
	}
}

void Compute(int& N, real* DI, real* AL, real* AU, int* IA, real* b)
{
	real *y = b, *x = b;
	// Решение Ly = b прямым обходом
	for (int i = 0; i < N; i++)
	{
		real sumL = 0;
		for (int k = 1; k < IA[i + 1] - IA[i] + 1; k++)
			sumL += AL[IA[i + 1] - k] * y[i - k];
		y[i] = (b[i] - sumL) / DI[i];
	}

	// Решение Ux = y обратным обходом
	for (int i = N - 1; i >= 0; i--)
	{
		real sumU = 0;
		int height = 1;

		for (int k = i + 1; k < N; k++)
		{
			if (IA[k + 1] - IA[k] >= height)
				sumU += x[k] * AU[IA[k + 1] - height];

			height++;
		}

		x[i] = y[i] - sumU;
	}
}

void Print(int& N, int& ALSize, real *DI, real *AL, real *AU, real *x)
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