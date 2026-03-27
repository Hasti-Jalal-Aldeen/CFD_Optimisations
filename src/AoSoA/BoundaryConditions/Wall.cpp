#include "AoSoA/BoundaryConditions.hpp"
using namespace AoSoA;

void Wall::update_cells(double* U[NVARS]) {
    int *CL, *CR;
    
    for (int fb = FB1; fb < FB2; fb++) {
        CL = domain.C[fb][0];
        CR = domain.C[fb][1];
        
        #pragma omp simd simdlen(VECLEN) safelen(VECLEN)
        for (int v = 0; v < VECLEN; v++) {
            U[0][CR[v]] = U[0][CL[v]];
        }
        for (int var = 1; var < NDIMS+1; var++) {
            #pragma omp simd simdlen(VECLEN) safelen(VECLEN)
            for (int v = 0; v < VECLEN; v++) {
                U[var][CR[v]] = -U[var][CL[v]];
            }
        }
        #pragma omp simd simdlen(VECLEN) safelen(VECLEN)
        for (int v = 0; v < VECLEN; v++) {
            U[NDIMS+1][CR[v]] = U[NDIMS+1][CL[v]];
        }
    }
}

