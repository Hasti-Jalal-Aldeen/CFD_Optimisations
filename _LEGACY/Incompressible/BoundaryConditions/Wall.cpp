#include "Incompressible/BoundaryConditions.hpp"

void Wall::update_cells() {
    for (int f = F1; f < F2; f++) {
        int CL = domain.C[f][0];
        int CR = domain.C[f][1];
        
        for (int dim = 0; dim < NDIMS; dim++) {
            U[dim][CR] = -U[dim][CL];
        }
        P[CR] = P[CL];
    }
}

void Wall::modify_mom_matrix(int N, double* A, double* rhs[NDIMS]) {
    for (int f = F1; f < F2; f++) {
        int CL = domain.C[f][0];
        int CR = domain.C[f][1];

        A[CR+CR*N] = 1.0;
        A[CR+CL*N] = 1.0;
    }

    for (int var = 0; var < nvars; var++) {
        for (int f = F1; f < F2; f++) {
            int CR = domain.C[f][1];
            rhs[var][CR] = 0.0;
        }
    }
}

void Wall::modify_pressure_matrix(int N, double* A, double* rhs) {
    modify_matrix_Neumann(F1, F2, N, A, domain);
    modify_rhs_Neumann_0(F1, F2, rhs, domain);
}
