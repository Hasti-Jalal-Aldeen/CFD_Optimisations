#include "Incompressible/BoundaryConditions.hpp"
#include <cmath>

void Symmetry::update_cells() {
    for (int f = F1; f < F2; f++) {
        int CL = domain.C[f][0];
        int CR = domain.C[f][1];
        
        double A[NDIMS];
        for (int dim = 0; dim < NDIMS; dim++) {
            A[dim] = domain.Af[f][dim];
        }
        double iA = 1.0/vabs(A);
        double n[NDIMS];
        for (int dim = 0; dim < NDIMS; dim++) {
            n[dim] = A[dim]*iA;
        }

        double ul[NDIMS];
        for (int dim = 0; dim < NDIMS; dim++) {
            ul[dim] = U[dim][CL];
        }
        
        double nor = dot(ul, n);
        for (int dim = 0; dim < NDIMS; dim++) {
            U[dim][CR] = ul[dim] - 2*nor*n[dim];
        }
        P[CR] = P[CL];
    }
}

void Symmetry::modify_mom_matrix(int N, double* A, double* rhs[NDIMS]) {
    // NEEDS CORRECTING
    modify_matrix_Dirichlet(F1, F2, N, A, domain);
    for (int var = 0; var < nvars; var++) {
        modify_rhs_Dirichlet(F1, F2, rhs[var], UBC[var], domain);
    }
}

void Symmetry::modify_pressure_matrix(int N, double* A, double* rhs) {
    modify_matrix_Neumann(F1, F2, N, A, domain);
    modify_rhs_Neumann_0(F1, F2, rhs, domain);
}
