#include "Base/GradientSchemes.hpp"
using namespace Base;

GradientScheme::GradientScheme(const Domain& domain)
: Nc(domain.Nc), Nf(domain.Nf), domain(domain)
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

    for (int f = 0; f < Nf; f++) {
        int CL = domain.C[f][0];
        int CR = domain.C[f][1];

        for (int var = 0; var < nvars; var++) {
            for (int dim = 0; dim < NDIMS; dim++) {
                double uf = u[var][CL]*domain.If[f] + u[var][CR]*(1-domain.If[f]);
            
                grad_u[var][dim*Nc+CL] += uf*domain.Af[f]*domain.nf[f][dim];
                grad_u[var][dim*Nc+CR] -= uf*domain.Af[f]*domain.nf[f][dim];
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
