#include "AoSoA/BoundaryConditions.hpp"
using namespace AoSoA;

void SupInlet::update_cells(double* U[NVARS]) {
    int *CL, *CR;

    for (int fb = FB1; fb < FB2; fb++) {
        CL = domain.C[fb][0];
        CR = domain.C[fb][1];

        for (int var = 0; var < nvars; var++) {
            #pragma omp simd simdlen(VECLEN) safelen(VECLEN)
            for (int v = 0; v < VECLEN; v++) {
                U[var][CR[v]] = UBC[var];
            }
        }
    }
}

void SubInlet::update_cells(double* U[NVARS]) {
    int *CL, *CR;
    double UL [NVARS][VECLEN] ALIGNED;
    double qql[VECLEN]        ALIGNED;
    double PL [VECLEN]        ALIGNED;

    for (int fb = FB1; fb < FB2; fb++) {
        CL = domain.C[fb][0];
        CR = domain.C[fb][1];

        for (int var = 0; var < nvars; var++) {
            #pragma omp simd simdlen(VECLEN) safelen(VECLEN)
            for (int v = 0; v < VECLEN; v++) {
                UL[var][v] = U[var][CL[v]];
            }
        }
        dotv(&UL[1], &UL[1], qql);
        #pragma omp simd simdlen(VECLEN) safelen(VECLEN)
        for (int v = 0; v < VECLEN; v++) {
            PL[v] = (gma-1)*(UL[NDIMS+1][v] - 0.5*qql[v]/UL[0][v]);
        }

        for (int var = 0; var < nvars; var++) {
            #pragma omp simd simdlen(VECLEN) safelen(VECLEN)
            for (int v = 0; v < VECLEN; v++) {
                U[var][CR[v]] = UBC[var];
            }
        }
        #pragma omp simd simdlen(VECLEN) safelen(VECLEN)
        for (int v = 0; v < VECLEN; v++) {
            U[NDIMS+1][CR[v]] = PL[v]/(gma-1) + 0.5*QBC[0]*dot(&QBC[1], &QBC[1]);
        }
    }
}
