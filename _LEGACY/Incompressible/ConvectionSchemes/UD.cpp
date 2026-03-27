#include "Incompressible/ConvectionSchemes.hpp"

void UD::discretise_momentum(double* A, double* b[NDIMS], double* F, double* G[NDIMS], const Domain& domain) const {
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
}
