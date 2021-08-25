#include "../stdafx.h"
#include "../Matrix.h"


namespace Solvers
{
	int LOS(SparseMatrix& A, std::vector<double>& x, std::vector<double>& b);
	int BCG(SparseMatrix& A, std::vector<double>& x, std::vector<double>& b);
}
