#include "Base/BoundaryConditions.hpp"
using namespace Base;

void SupInlet::update_cells(double* U[NVARS]) {
    for (int f = F1; f < F2; f++) {
        int CL = domain.C[f][0];
        int CR = domain.C[f][1];

        for (int var = 0; var < nvars; var++) {
            U[var][CR] = UBC[var];
        }
    }
}

void SubInlet::update_cells(double* U[NVARS]) {
    double UL[nvars];
    for (int f = F1; f < F2; f++) {
        int CL = domain.C[f][0];
        int CR = domain.C[f][1];

        for (int var = 0; var < nvars; var++) {
            UL[var] = U[var][CL];
        }
        double PL = (gma-1)*(UL[NDIMS+1] - 0.5*dot(&UL[1], &UL[1])/UL[0]);

        for (int var = 0; var < nvars; var++) {
            U[var][CR] = UBC[var];
        }
        U[NDIMS+1][CR] = PL/(gma-1) + 0.5*QBC[0]*dot(&QBC[1], &QBC[1]);
    }
}
