#include "Incompressible/ConvectionSchemes.hpp"

void LUD::discretise_momentum(double* A, double* b[NDIMS], double* F, double* G[NDIMS], const Domain& domain) const {
    int Nc = domain.Nc;
    int Nf = domain.Nf;

    for (int f = 0; f < Nf; f++) {
        int CL = domain.C[f][0];
        int CR = domain.C[f][1];

        double Ff = F[f];
        int CU = (Ff >= 0) ? CL : CR; // upwind cell

        A[CL+CU*Nc] += Ff*domain.iv[CL];
        A[CR+CU*Nc] -= Ff*domain.iv[CR];
    }

    for (int var = 0; var < NDIMS; var++) {
        for (int f = 0; f < Nf; f++) {
            double Ff = F[f];
            int LR = (Ff < 0);
            
            int CL = domain.C[f][0];
            int CR = domain.C[f][1];
            int CU = (Ff >= 0) ? CL : CR; // upwind cell
    
            double r[NDIMS], GU[NDIMS];
            for (int dim = 0; dim < NDIMS; dim++) {
                GU[dim] = G[var][dim*Nc+CU];
                r[dim]  = domain.r[f][LR][dim];
            }

            double r_dot_G = dot(r, GU);

            b[var][CL] -= Ff*r_dot_G*domain.iv[CL];
            b[var][CR] += Ff*r_dot_G*domain.iv[CR];
        }
    }
}