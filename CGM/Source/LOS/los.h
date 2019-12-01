#include "../Matrix/matrix.h"
#include <cmath>

int LOS(Matrix& A, double* x, double* f, AuxVectors& aux, int maxiter, double eps, double& lastdiff);
int LOS_diag1(Matrix& A, double* x, double* f, AuxVectors& aux, int maxiter, double eps, double& lastdiff);
int LOS_diag2(Matrix& A, double* x, double* f, AuxVectors& aux, int maxiter, double eps, double& lastdiff);
int LOS_diag3(Matrix& A, double* x, double* f, AuxVectors& aux, int maxiter, double eps, double& lastdiff);
int LOS_LU(Matrix& A, double* x, double* f, Matrix& LU, AuxVectors& aux, int maxiter, double eps, double& lastdiff);