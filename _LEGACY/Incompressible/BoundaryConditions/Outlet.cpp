#include "Incompressible/BoundaryConditions.hpp"

void Outlet::update_cells() {
    for (int f = F1; f < F2; f++) {
        int CL = domain.C[f][0];
        int CR = domain.C[f][1];

        for (int var = 0; var < nvars; var++) {
            U[var][CR] = U[var][CL];
        }
        P[CR] = P_BC;
    }
}

void Outlet::modify_mom_matrix(int N, double* A, double* rhs[NDIMS]) {
    modify_matrix_Neumann(F1, F2, N, A, domain);
    for (int var = 0; var < nvars; var++) {
        modify_rhs_Neumann_0(F1, F2, rhs[var], domain);
    }
}

void Outlet::modify_pressure_matrix(int N, double* A, double* rhs) {
    modify_matrix_Dirichlet(F1, F2, N, A, domain);
    modify_rhs_Dirichlet(F1, F2, rhs, P_BC, domain);
}
