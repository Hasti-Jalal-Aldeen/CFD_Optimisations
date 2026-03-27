#include "SoA/GradientSchemes.hpp"
using namespace SoA;

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
        int CL = domain.C[0][f];
        int CR = domain.C[1][f];
        double A  = domain.Af[f];

        for (int var = 0; var < nvars; var++) {
            for (int dim = 0; dim < NDIMS; dim++) {
                double uf = u[var][CL]*domain.If[f] + u[var][CR]*(1-domain.If[f]);

                grad_u[var][dim*Nc+CL] += uf*A*domain.nf[dim][f];
                grad_u[var][dim*Nc+CR] -= uf*A*domain.nf[dim][f];
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