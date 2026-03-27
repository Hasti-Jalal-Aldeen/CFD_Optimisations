#include "IterativeSolvers.hpp"

IterativeSolver::IterativeSolver(const int maxniters, const int N, const double tol) :
maxniters(maxniters), N(N), tol(tol) 
{
    return;
}

int IterativeSolver::solve(const double* A, double* x, const double* b) const {
    return -1;
}