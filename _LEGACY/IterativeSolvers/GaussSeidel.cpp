#include "IterativeSolvers.hpp"
#include <math.h>

int GaussSeidel::solve(const double* A, double* x, const double* b) const {
    int iter;
    for (iter = 0; iter < maxniters; iter++) {
        if (iter % 20 == 0) {
            // residual
            double res = 0.0;
            for (int i = 0; i < N; i++) {
                double row = b[i];
                for (int j = 0; j < N; j++) {
                    row -= A[i+j*N]*x[j];
                }
                res += row*row;
            }
            res = sqrt(res);
            if (res < tol) return iter;
        }

        // gauss seidel
        // forward pass
        for (int i = 0; i < N; i++) {
            double term1 = 0.0;
            for (int j = 0; j < i; j++) {
                term1 += A[i+j*N]*x[j];
            }

            double term2 = 0.0;
            for (int j = i+1; j < N; j++) {
                term2 += A[i+j*N]*x[j];
            }

            x[i] = (b[i] - term1 - term2)/A[i+i*N];
        }
        // backward pass
        for (int i = N-1; i >= 0; i--) {
            double term1 = 0.0;
            for (int j = 0; j < i; j++) {
                term1 += A[i+j*N]*x[j];
            }

            double term2 = 0.0;
            for (int j = i+1; j < N; j++) {
                term2 += A[i+j*N]*x[j];
            }

            x[i] = (b[i] - term1 - term2)/A[i+i*N];
        }
    }
    return iter;
}