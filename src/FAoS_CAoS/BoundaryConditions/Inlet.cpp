#include "FAoS_CAoS/BoundaryConditions.hpp"
using namespace FAoS_CAoS;

void SupInlet::update_cells(double U[][NVARS]) {
    for (int f = F1; f < F2; f++) {
        int CL = domain.C[f][0];
        int CR = domain.C[f][1];

        for (int var = 0; var < nvars; var++) {
            U[CR][var] = UBC[var];
        }
    }
}

void SubInlet::update_cells(double U[][NVARS]) {
    double UL[nvars];
    for (int f = F1; f < F2; f++) {
        int CL = domain.C[f][0];
        int CR = domain.C[f][1];

        double PL = (gma-1)*(U[CL][NDIMS+1] - 0.5*dot(&U[CL][1], &U[CL][1])/U[CL][0]);

        for (int var = 0; var < NVARS; var++) {
            U[CR][var] = UBC[var];
        }
        U[CR][NDIMS+1] = PL/(gma-1) + 0.5*QBC[0]*dot(&QBC[1], &QBC[1]);
    }
}
