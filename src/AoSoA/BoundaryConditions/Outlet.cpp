#include "AoSoA/BoundaryConditions.hpp"
using namespace AoSoA;

void SupOutlet::update_cells(double* U[NVARS]) {
    int *CL, *CR;

    for (int fb = FB1; fb < FB2; fb++) {
        CL = domain.C[fb][0];
        CR = domain.C[fb][1];

        for (int var = 0; var < nvars; var++) {
            #pragma omp simd simdlen(VECLEN) safelen(VECLEN)
            for (int v = 0; v < VECLEN; v++) {
                U[var][CR[v]] = U[var][CL[v]];
            }
        }
    }
}

void SubOutlet::update_cells(double* U[NVARS]) {
    int *CL, *CR;
    double rl[VECLEN]         ALIGNED;
    double rul[NDIMS][VECLEN] ALIGNED;
    double qql[VECLEN]        ALIGNED;

    for (int fb = FB1; fb < FB2; fb++) {
        CL = domain.C[fb][0];
        CR = domain.C[fb][1];

        #pragma omp simd simdlen(VECLEN) safelen(VECLEN)
        for (int v = 0; v < VECLEN; v++) {
            rl[v] = U[0][CL[v]];
        }
        for (int dim = 0; dim < NDIMS; dim++) {
            #pragma omp simd simdlen(VECLEN) safelen(VECLEN)
            for (int v = 0; v < VECLEN; v++) {
                rul[dim][v] = U[dim+1][CL[v]];
            }
        }

        for (int var = 0; var < NDIMS; var++) {
            #pragma omp simd simdlen(VECLEN) safelen(VECLEN)
            for (int v = 0; v < VECLEN; v++) {
                U[var][CR[v]] = rul[var-1][v];
            }
        }

        dotv(rul, rul, qql);
        #pragma omp simd simdlen(VECLEN) safelen(VECLEN)
        for (int v = 0; v < VECLEN; v++) {
            U[NDIMS+1][CR[v]] = QBC[NDIMS+1]/(gma-1) + 0.5*qql[v]/rl[v];
        }
    }
}
