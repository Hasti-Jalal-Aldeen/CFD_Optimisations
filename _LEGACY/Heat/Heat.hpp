#pragma once
#include "Global.hpp"
#include "Domain.hpp"
#include "IterativeSolvers.hpp"
#include "Writer.hpp"

class BC {
    public:
    BC();
    void init(int F1, int F2, double UBC, int C[][2], double Af[][NDIMS], double* U);
    virtual void modify_matrix(int N, double* A, double* rhs);
    virtual void update_cells();

    int F1, F2;
    double UBC;

    int (*C)[2];
    double (*Af)[NDIMS];

    double* U;
};

class Dirichlet : public BC {
    using BC::BC;
    void modify_matrix(int N, double* A, double* rhs) override;
    void update_cells() override;
};

class Neumann : public BC {
    using BC::BC;
    void modify_matrix(int N, double* A, double* rhs) override;
    void update_cells() override;
};

class Heat {
    public:
    Heat(const Domain& domain, const IterativeSolver& solver);
    ~Heat();
    void initialise(double TIC);
    void gather_BCs(const int BCmap[], const double TBC[]);
    void run_solver();

    private:
    Writer* writer;

    const int Nc, Nf;
    const Domain& domain;
    const IterativeSolver& solver;

    // BCs
    const int NBC;
    BC** BCs;
    void update_BC_cells();
    void update_BC_matrix();

    // cells
    double* T;

    // Heat Equation
    double* A;
    double* rhs;
};