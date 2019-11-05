#include "../Matrix/matrix.h"
#include "../Memory/memory.h"
#include <cmath>

int CGM(Matrix& A, double *x, double *f, Memory& cache, int maxiter, double eps);
int CGM_diag(Matrix& A, double *x, double *f, Memory& cache, int maxiter, double eps);