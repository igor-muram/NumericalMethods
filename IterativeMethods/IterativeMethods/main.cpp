#include "Matrix.h"

int main()
{
	int N = 0, di_num = 0, m = 0;
	InputSize(N, di_num, m);

	real *b = new real[N];
	real *res = new real[N];

	real **A = new real* [N];
	for (int i = 0; i < N; i++)
		A[i] = new real[di_num];

	for (int i = 0; i < N; i++)
		res[i] = 0;

	Input(N, di_num, A, b);
	Multiply(N, di_num, A, m, b, res);
	Output(N, res);

	system("pause");
	return 0;
}