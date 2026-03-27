#include "Base/BoundaryConditions.hpp"
using namespace Base;

void SupOutlet::update_cells(double* U[NVARS]) {
    for (int f = F1; f < F2; f++) {
        int CL = domain.C[f][0];
        int CR = domain.C[f][1];

        for (int var = 0; var < nvars; var++) {
            U[var][CR] = U[var][CL];
        }
    }
}

void SubOutlet::update_cells(double* U[NVARS]) {
    for (int f = F1; f < F2; f++) {
        int CL = domain.C[f][0];
        int CR = domain.C[f][1];

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
