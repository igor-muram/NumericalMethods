#include "../Matrix/matrix.h"
#include "../Memory/memory.h"
#include <cmath>

int LOS(Matrix& A, double *x, double *f, Memory& cache, int maxiter, double eps);
int LOS_diag(Matrix& A, double *x, double *f, Memory& cache, int maxiter, double eps);
