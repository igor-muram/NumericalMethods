#include "Matrix.h"

int main()
{
	int N = 0, ALSize = 0;
	InputSize(N, ALSize);

	real *DI = new real[N];
	real *AL = new real[ALSize];
	real *AU = new real[ALSize];
	real *b = new real[N];
	real *x = new real[N];
	int *IA = new int[N + 1];

	real **A = new real*[N];
	for (int i = 0; i < N; i++)
		A[i] = new real[N];
	
	Input(N, ALSize, DI, AL, AU, b, IA);
	BuildLU(N, DI, AL, AU, IA);
	Compute(N, DI, AL, AU, IA, b);
	Output(N, ALSize, DI, AL, AU, b);

	system("pause");
	return 0;
}