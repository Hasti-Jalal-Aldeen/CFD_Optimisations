#include "AoSoA/GradientSchemes.hpp"
using namespace AoSoA;

GradientScheme::GradientScheme(const Domain& domain)
: Nc(domain.Nc), Nfb(domain.Nfb), domain(domain)
{
    return;
}

void GradientScheme::compute_gradient(int nvars, const double* const u[NVARS], double* grad_u[NVARS]) const {
    return;
}

void GaussGreen::compute_gradient(int nvars, const double* const u[NVARS], double* grad_u[NVARS]) const {
    for (int var = 0; var < nvars; var++) {
        setarr(0, NDIMS*Nc, 0.0, grad_u[var]);
    }

    int *CL, *CR;
    double uf[VECLEN] ALIGNED;

    for (int fb = 0; fb < Nfb; fb++) {
        CL = domain.C[fb][0];
        CR = domain.C[fb][1];

        for (int var = 0; var < nvars; var++) {
            for (int dim = 0; dim < NDIMS; dim++) {
                for (int v = 0; v < VECLEN; v++) {
                    uf[v] = u[var][CL[v]]*domain.If[fb][v] + u[var][CR[v]]*(1-domain.If[fb][v]);
                }
            
                for (int v = 0; v < VECLEN; v++) {
                    grad_u[var][dim*Nc+CL[v]] += uf[v]*domain.Af[fb][v]*domain.nf[fb][dim][v];
                    grad_u[var][dim*Nc+CR[v]] -= uf[v]*domain.Af[fb][v]*domain.nf[fb][dim][v];
                }
            }
        }
    }
    
    for (int var = 0; var < nvars; var++) {
        for (int dim = 0; dim < NDIMS; dim++) {
            for (int i = 0; i < Nc; i++) {
                grad_u[var][dim*Nc+i] *= domain.iv[i];
            }
        }
    }
}
