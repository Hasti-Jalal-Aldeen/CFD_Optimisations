#include "VEC/GradientSchemes.hpp"
using namespace VEC;

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

    int     *CL, *CR;
    double  *Af, *If;
    double  nf[NDIMS][VECLEN] ALIGNED;
    double  uf[VECLEN]        ALIGNED;

    for (int fb = 0; fb < Nf; fb+=VECLEN) {
        CL = &domain.C[0][fb];
        CR = &domain.C[1][fb];

        for (int dim = 0; dim < NDIMS; dim++) {
            for (int v = 0; v < VECLEN; v++) {
                nf[dim][v] = domain.nf[dim][fb+v];
            }
        }

        Af = &domain.Af[fb];
        If = &domain.If[fb];

        for (int var = 0; var < nvars; var++) {
            for (int v = 0; v < VECLEN; v++) {
                uf[v] = u[var][CL[v]]*If[v] + u[var][CR[v]]*(1-If[v]);
            }

            for (int dim = 0; dim < NDIMS; dim++) {
                for (int v = 0; v < VECLEN; v++) {
                    double val = uf[v]*Af[v]*nf[dim][v];
                    grad_u[var][dim*Nc+CL[v]] += val;
                    grad_u[var][dim*Nc+CR[v]] -= val;
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