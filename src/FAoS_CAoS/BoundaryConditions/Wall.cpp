#include "FAoS_CAoS/BoundaryConditions.hpp"
using namespace FAoS_CAoS;

void Wall::update_cells(double U[][NVARS]) {
    for (int f = F1; f < F2; f++) {
        int CL = domain.C[f][0];
        int CR = domain.C[f][1];
        
        U[CR][0] = U[CL][0];
        for (int var = 1; var < NDIMS+1; var++) {
            U[CR][var] = -U[CL][var];
        }
        U[CR][NDIMS+1] = U[CL][NDIMS+1];
    }
}

