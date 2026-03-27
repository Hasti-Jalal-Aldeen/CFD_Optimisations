#include "SoA/BoundaryConditions.hpp"
using namespace SoA;

BoundaryCondition::BoundaryCondition(const Domain& domain) : domain(domain) {
    Nc = domain.Nc;
}

void BoundaryCondition::init(int nvars, int F1, int F2, double gma, const double QBC[]) {
    this->nvars = nvars;
    this->F1 = F1;
    this->F2 = F2;
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
    for (int var = 0; var < nvars; var++) {
        for (int dim = 0; dim < NDIMS; dim++) {
            for (int f = F1; f < F2; f++) {
                int CL = domain.C[0][f];
                int CR = domain.C[1][f];
                
                G[var][dim*Nc+CR] = -G[var][dim*Nc+CL];
            }
        }
    }
}

void BoundaryCondition::update_limiter(double* phi[NVARS]) {
    const int Nc = domain.Nc;
    for (int var = 0; var < nvars; var++) {
        for (int f = F1; f < F2; f++) {
            int CL = domain.C[0][f];
            int CR = domain.C[1][f];
            
            phi[var][CR] = phi[var][CL];
        }
    }
}
