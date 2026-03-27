#include "IterativeSolvers.hpp"
#include "cblas.h"

ConjugateGradient::ConjugateGradient(const int maxniters, const int N, const double tol) : IterativeSolver(maxniters, N, tol) {
    r  = new double[N];
    p  = new double[N];
    Ap = new double[N];
}

ConjugateGradient::~ConjugateGradient() {
    delete[] r;
    delete[] p;
    delete[] Ap;
}

int ConjugateGradient::solve(const double* A, double* x, const double* b) const {
    cblas_dcopy(N, b, 1, r, 1);
    cblas_dgemv(CblasColMajor, CblasNoTrans, N, N, -1.0, A, N, x, 1, 1.0, r, 1);
    cblas_dcopy(N, r, 1, p, 1);

    double rkp1_dot_rkp1 = cblas_ddot(N, r, 1, r, 1);
    int iter;
    for (iter = 0; iter < maxniters; iter++) {
        double rk_dot_rk = rkp1_dot_rkp1;
        cblas_dgemv(CblasColMajor, CblasNoTrans, N, N, 1.0, A, N, p, 1, 0.0, Ap, 1);
        double p_dot_Ap = cblas_ddot(N, p, 1, Ap, 1);
        double alpha = rk_dot_rk/p_dot_Ap;
        cblas_daxpy(N,  alpha,  p, 1, x, 1);
        cblas_daxpy(N, -alpha, Ap, 1, r, 1);
        double nrm_r = cblas_dnrm2(N, r, 1);
        if (nrm_r < tol) return iter;
        rkp1_dot_rkp1 = cblas_ddot(N, r, 1, r, 1);
        double beta = rkp1_dot_rkp1/rk_dot_rk;
        cblas_dscal(N, beta, p, 1);
        cblas_daxpy(N, 1.0, r, 1, p, 1);
    }
    return iter;
}
