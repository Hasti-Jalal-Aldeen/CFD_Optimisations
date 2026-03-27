#include "CAoS/BoundaryConditions.hpp"
using namespace CAoS;

void SupOutlet::update_cells(double U[][NVARS]) {
    for (int f = F1; f < F2; f++) {
        int CL = domain.C[0][f];
        int CR = domain.C[1][f];

        for (int var = 0; var < NVARS; var++) {
            U[CR][var] = U[CL][var];
        }
    }
}

void SubOutlet::update_cells(double U[][NVARS]) {
    for (int f = F1; f < F2; f++) {
        int CL = domain.C[0][f];
        int CR = domain.C[1][f];

        for (int var = 0; var < NVARS; var++) {
            U[CR][var] = U[CL][var];
        }

        U[CR][NDIMS+1] = QBC[NDIMS+1]/(gma-1) + 0.5*dot(&U[CL][1],&U[CL][1])/U[CL][0];
    }
}
