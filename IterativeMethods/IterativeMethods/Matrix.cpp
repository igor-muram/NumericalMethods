#include "Matrix.h"

void InputSize(int& N, int& di_num, int& m)
{
	ifstream in("input/size.txt", ios::in);
	in >> N  >> di_num >> m;
	in.close();
}

void Input(int& N, int& di_num, real **A, real *b)
{
	ifstream in("input/matrix.txt");
	for (int i = 0; i < N; i++)
		for (int j = 0; j < di_num; j++)
			in >> A[i][j];
	in.close();

	in.open("input/b.txt");
	for (int i = 0; i < N; i++)
		in >> b[i];
	in.close();
}

void Multiply(int& N, int& di_num, real **A, int& m, real *b, real *res)
{
	int i_di = (di_num - 1) / 2;
	int i_m = m + 1;

	for (int i = 0; i < N; i++)
		res[i] += A[i][i_di] * b[i];

	for (int i = 1; i < N; i++)
	{
		res[i] += A[i][i_di - 1] * b[i - 1];
		res[i - 1] += A[i][i_di + 1] * b[i];
	}

	for (int i = i_m; i < N; i++)
	{
		res[i] += A[i][i_di - 2] * b[i - i_m];
		res[i - i_m] += A[i][i_di + 2] * b[i];
	}
}

void Output(int& N, real *res)
{
	cout << "Result: " << endl;
	for (int i = 0; i < N; i++)
		cout << res[i] << ' ';
	cout << endl << endl;
}