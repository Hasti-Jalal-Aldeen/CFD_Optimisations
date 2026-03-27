#include "FAoS_CAoS/GradientSchemes.hpp"
using namespace FAoS_CAoS;

GradientScheme::GradientScheme(const Domain& domain)
: Nc(domain.Nc), Nf(domain.Nf), domain(domain)
{
    return;
}

void GradientScheme::compute_gradient(int nvars, const double u[][NVARS], double grad_u[][NVARS][NDIMS]) const {
    return;
}

void GaussGreen::compute_gradient(int nvars, const double u[][NVARS], double grad_u[][NVARS][NDIMS]) const {
    setarr(0, nvars*NDIMS*Nc, 0.0, grad_u[0][0]);

    for (int f = 0; f < Nf; f++) {
        int CL = domain.C[f][0];
        int CR = domain.C[f][1];
        double A = domain.Af[f];

        for (int var = 0; var < NVARS; var++) {
            for (int dim = 0; dim < NDIMS; dim++) {
                double uf = u[CL][var]*domain.If[f] + u[CR][var]*(1-domain.If[f]);

                grad_u[CL][var][dim] += uf*A*domain.nf[f][dim];
                grad_u[CR][var][dim] -= uf*A*domain.nf[f][dim];
            }
        }
    }
    
    for (int i = 0; i < Nc; i++) {
        for (int var = 0; var < NVARS; var++) {
            for (int dim = 0; dim < NDIMS; dim++) {
                grad_u[i][var][dim] *= domain.iv[i];
            }
        }
    }
}