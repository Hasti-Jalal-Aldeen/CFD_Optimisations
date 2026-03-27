#include "Base/BoundaryConditions.hpp"
using namespace Base;

void Wall::update_cells(double* U[NVARS]) {
    for (int f = F1; f < F2; f++) {
        int CL = domain.C[f][0];
        int CR = domain.C[f][1];
        
        U[0][CR] = U[0][CL];
        for (int var = 1; var < NDIMS+1; var++) {
            U[var][CR] = -U[var][CL];
        }
        U[NDIMS+1][CR] = U[NDIMS+1][CL];
    }
}

