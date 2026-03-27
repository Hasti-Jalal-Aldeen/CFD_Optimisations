#include "AoSoA/BoundaryConditions.hpp"
#include <cmath>
using namespace AoSoA;

void Symmetry::update_cells(double* U[NVARS]) {
    int *CL, *CR;
    double nf[NDIMS][VECLEN] ALIGNED;
    double ul[NDIMS][VECLEN] ALIGNED;
    double nor[VECLEN]       ALIGNED;

    for (int fb = FB1; fb < FB2; fb++) {
        CL = domain.C[fb][0];
        CR = domain.C[fb][1];
        
        for (int dim = 0; dim < NDIMS; dim++) {
            #pragma omp simd simdlen(VECLEN) safelen(VECLEN)
            for (int v = 0; v < VECLEN; v++) {
                nf[dim][v] = domain.nf[fb][dim][v];
            }
        }

        for (int dim = 0; dim < NDIMS; dim++) {
            #pragma omp simd simdlen(VECLEN) safelen(VECLEN)
            for (int v = 0; v < VECLEN; v++) {
                ul[dim][v] = U[dim+1][CL[v]];
            }
        }
        dotv(ul, nf, nor);

        #pragma omp simd simdlen(VECLEN) safelen(VECLEN)
        for (int v = 0; v < VECLEN; v++) {
            U[0][CR[v]] = U[0][CL[v]];
        }
        for (int dim = 0; dim < NDIMS; dim++) {
            #pragma omp simd simdlen(VECLEN) safelen(VECLEN)
            for (int v = 0; v < VECLEN; v++) {
                U[dim+1][CR[v]] = ul[dim][v] - 2*nor[v]*nf[dim][v];
            }
        }
        #pragma omp simd simdlen(VECLEN) safelen(VECLEN)
        for (int v = 0; v < VECLEN; v++) {
            U[NDIMS+1][CR[v]] = U[NDIMS+1][CL[v]];
        }
    }
}
