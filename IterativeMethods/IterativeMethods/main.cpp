#include "Matrix.h"

int main()
{
	Matrix mat = { };
	ReadMatrix(mat);

	real* xstar = new real[mat.N];
	for (int i = 0; i < mat.N; i++)
		xstar[i] = i + 1;

	real* f = new real[mat.N];
	for (int i = 0; i < mat.N; i++)
		f[i] = Multiply(mat, xstar, i);
	
	real* x0 = new real[mat.N];
	ReadX0(x0, mat.N);

	real* x = new real[mat.N];
	real* r = new real[mat.N];

	for (int i = 90; i <= 100; i++)
	{
		ReadX0(x0, mat.N);
		for (int i = 0; i < mat.N; i++)
			x[i] = 0;

		real w = 0.01 * i;
		int J = Jacobi(mat, f, x, x0, r, w);

		real sum1 = 0, sum2 = Norm(x, mat.N), sum3 = 0, sum4 = Norm(f, mat.N);
		for (int i = 0; i < mat.N; i++)
		{
			sum1 += (xstar[i] - x[i]) * (xstar[i] - x[i]);
			real temp = Multiply(mat, x, i);
			sum3 += (temp - f[i]) * (temp - f[i]);
		}

		sum1 = sqrt(sum1);
		sum3 = sqrt(sum3);
		real cond = (sum1 * sum4) / (sum2 * sum3);

		printf_s("%.2f\t\t%.15f\t%.15f\t%d\t%.6f\n", w, x[0], fabs(xstar[0] - x[0]), J, cond);

		for (int i = 1; i < mat.N; i++)
			printf_s("\t\t%.15f\t%.15f\n", x[i], fabs(xstar[i] - x[i]));

		printf_s("\n\n");
	}

	system("pause");
	return 0;
}