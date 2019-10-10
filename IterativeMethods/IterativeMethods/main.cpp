#include "Matrix.h"

int main()
{
	Matrix mat = { };
	ReadMatrix(mat);

	real *f = new real[mat.N];
	real *x0 = new real[mat.N];
	real *x = new real[mat.N];
	real *r = new real[mat.N];

	ReadF(f, mat.N);
	ReadX0(x0, mat.N);

	for (int i = 0; i < mat.N; i++)
		x[i] = x0[i];

	int J = Jacobi(mat, f, x0, x, r, 0.999);
	ReadX0(x0, mat.N);
	int Z = Zeidel(mat, f, x0, r, 1);

	system("pause");
	return 0;
}