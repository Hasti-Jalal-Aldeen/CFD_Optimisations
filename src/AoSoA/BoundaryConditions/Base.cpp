#include "AoSoA/BoundaryConditions.hpp"
using namespace AoSoA;

BoundaryCondition::BoundaryCondition(const Domain& domain) : domain(domain) {
    Nc = domain.Nc;
}

void BoundaryCondition::init(int nvars, int FB1, int FB2, double gma, const double QBC[]) {
    this->nvars = nvars;
    this->FB1 = FB1;
    this->FB2 = FB2;
    this->gma = gma;

    for (int var = 0; var < nvars; var++) {
        this->QBC[var] = QBC[var];
    }
    
    UBC[0] = QBC[0];
    for (int var = 1; var < nvars; var++) {
        UBC[var] = QBC[0]*QBC[var];
    }
    UBC[NDIMS+1] = QBC[NDIMS+1]/(gma-1) + 0.5*QBC[0]*dot(&QBC[1], &QBC[1]);
}

void BoundaryCondition::update_cells(double* U[NVARS]) {
    return;
}

void BoundaryCondition::update_grads(double* G[NVARS]) {
    const int Nc = domain.Nc;
    int *CL, *CR;

    for (int var = 0; var < nvars; var++) {
        for (int dim = 0; dim < NDIMS; dim++) {
            for (int fb = FB1; fb < FB2; fb++) {
                CL = domain.C[fb][0];
                CR = domain.C[fb][1];
                
                #pragma omp simd simdlen(VECLEN) safelen(VECLEN)
                for (int v = 0; v < VECLEN; v++) {
                    G[var][dim*Nc+CR[v]] = -G[var][dim*Nc+CL[v]];
                }
            }
        }
    }
}

void BoundaryCondition::update_limiter(double* phi[NVARS]) {
    const int Nc = domain.Nc;
    int *CL, *CR;

    for (int var = 0; var < nvars; var++) {
        for (int fb = FB1; fb < FB2; fb++) {
            CL = domain.C[fb][0];
            CR = domain.C[fb][1];
            
            #pragma omp simd simdlen(VECLEN) safelen(VECLEN)
            for (int v = 0; v < VECLEN; v++) {
                phi[var][CR[v]] = phi[var][CL[v]];
            }
        }
    }
}
