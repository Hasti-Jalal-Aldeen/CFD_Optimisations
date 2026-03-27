#include "VEC/BoundaryConditions.hpp"
using namespace VEC;

void SupOutlet::update_cells(double* U[NVARS]) {
    for (int f = F1; f < F2; f++) {
        int CL = domain.C[0][f];
        int CR = domain.C[1][f];

        for (int var = 0; var < nvars; var++) {
            U[var][CR] = U[var][CL];
        }
    }
}

void SubOutlet::update_cells(double* U[NVARS]) {
    for (int f = F1; f < F2; f++) {
        int CL = domain.C[0][f];
        int CR = domain.C[1][f];

        double rl = U[0][CL];
        double rul[NDIMS];
        for (int dim = 0; dim < NDIMS; dim++) {
            rul[dim] = U[dim+1][CL];
        }

        for (int var = 0; var < NDIMS; var++) {
            U[var][CR] = rul[var-1];
        }

        U[NDIMS+1][CR] = QBC[NDIMS+1]/(gma-1) + 0.5*dot(rul,rul)/rl;
    }
}
