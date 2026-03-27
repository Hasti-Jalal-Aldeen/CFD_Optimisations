#include "SoA/BoundaryConditions.hpp"
using namespace SoA;

void Wall::update_cells(double* U[NVARS]) {
    for (int f = F1; f < F2; f++) {
        int CL = domain.C[0][f];
        int CR = domain.C[1][f];
        
        U[0][CR] = U[0][CL];
        for (int var = 1; var < NDIMS+1; var++) {
            U[var][CR] = -U[var][CL];
        }
        U[NDIMS+1][CR] = U[NDIMS+1][CL];
    }
}

