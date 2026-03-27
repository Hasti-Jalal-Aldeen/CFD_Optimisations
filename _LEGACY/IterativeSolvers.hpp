#pragma once

class IterativeSolver {
    public:
    IterativeSolver(const int maxniters, const int N, const double tol);
    virtual int solve(const double* A, double* x, const double* b) const;

    const int maxniters;
    const int N;
    const double tol;
};

class GaussSeidel : public IterativeSolver {
    public:
    using IterativeSolver::IterativeSolver;
    int solve(const double* A, double* x, const double* b) const override;
};

class ConjugateGradient : public IterativeSolver {
    public:
    ConjugateGradient(const int maxniters, const int N, const double tol);
    ~ConjugateGradient();
    int solve(const double* A, double* x, const double* b) const override;

    private:
    double* r;
    double* p;
    double* Ap;
};

