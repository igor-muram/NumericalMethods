#include "../Matrix/matrix.h"
#include <cmath>

int CGM(Matrix& A, double* x, double* f, AuxVectors& aux, int maxiter, double eps, double& lastdiff);
int CGM_diag1(Matrix& A, double* x, double* f, AuxVectors& aux, int maxiter, double eps, double& lastdiff);
int CGM_diag2(Matrix& A, double* x, double* f, AuxVectors& aux, int maxiter, double eps, double& lastdiff);
int CGM_diag3(Matrix& A, double* x, double* f, AuxVectors& aux, int maxiter, double eps, double& lastdiff);
int CGM_LU(Matrix& A, double* x, double* f, Matrix& LU, AuxVectors& aux, int maxiter, double eps, double& lastdiff);