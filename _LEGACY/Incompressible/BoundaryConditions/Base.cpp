#include "Incompressible/BoundaryConditions.hpp"

void modify_matrix_Dirichlet(int F1, int F2, int N, double* A, const Domain& domain) {
    for (int f = F1; f < F2; f++) {
        int CL = domain.C[f][0];
        int CR = domain.C[f][1];

        A[CR+CR*N] = 1.0;
        A[CR+CL*N] = 0.0;
    }
}
void modify_rhs_Dirichlet(int F1, int F2, double* rhs, double UBC, const Domain& domain) {
    for (int f = F1; f < F2; f++) {
        int CR = domain.C[f][1];
        rhs[CR] = UBC;
    }
}

void modify_matrix_Neumann(int F1, int F2, int N, double* A, const Domain& domain) {
    for (int f = F1; f < F2; f++) {
        int CL = domain.C[f][0];
        int CR = domain.C[f][1];

        A[CR+CR*N] =  1.0;
        A[CR+CL*N] = -1.0;
    }
}
void modify_rhs_Neumann_0(int F1, int F2, double* rhs, const Domain& domain) {
    for (int f = F1; f < F2; f++) {
        int CR = domain.C[f][1];
        rhs[CR] = 0.0;
    }
}
void modify_rhs_Neumann  (int F1, int F2, double* rhs, double du_dn, const Domain& domain) {
    for (int f = F1; f < F2; f++) {
        int CL = domain.C[f][0];
        int CR = domain.C[f][1];

        double d[NDIMS];
        for (int dim = 0; dim < NDIMS; dim++) {
            d[dim] = domain.r[f][0][dim] - domain.r[f][1][dim];
        }

        rhs[CR] = du_dn * vabs(d);
    }
}

BoundaryCondition::BoundaryCondition(const Domain& domain) : domain(domain) {
    return;
}

void BoundaryCondition::init(int nvars, int F1, int F2, const double UBC[], double P_BC, double* U[NDIMS], double* P) {
    this->nvars = nvars;
    this->F1 = F1;
    this->F2 = F2;

    this->P = P;
    this->P_BC = P_BC;

    for (int var = 0; var < nvars; var++) {
        this->U[var] = U[var];
        this->UBC[var] = UBC[var];
    }
}

void BoundaryCondition::modify_mom_matrix(int N, double* A, double* rhs[NDIMS]) {
    return;
}

void BoundaryCondition::modify_pressure_matrix(int N, double* A, double* rhs) {
    return;
}

void BoundaryCondition::update_cells() {
    return;
}

/*
void BoundaryCondition::modify_mom_matrix(int N, double* A, double* rhs[NDIMS]) {
    for (int f = F1; f < F2; f++) {
        int CL = C[f][0];
        int CR = C[f][1];

        A[CR+CL*N] = 0.0;
        A[CR+CR*N] = 1.0;
    }

    for (int var = 0; var < nvars; var++) {
        for (int f = F1; f < F2; f++) {
            int CR = C[f][1];
            b[var][CR] = U[var][CR];
        }
    }
}

void BoundaryCondition::modify_pressure_matrix(int N, double* A, double* rhs) {
    for (int f = F1; f < F2; f++) {
        int CL = C[f][0];
        int CR = C[f][1];

        A[CR+CR*N] =  1.0;
        A[CR+CL*N] = -1.0;

        rhs[CR] = 0.0;
    }
}
*/