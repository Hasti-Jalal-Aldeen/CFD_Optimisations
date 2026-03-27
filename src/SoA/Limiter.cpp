#include "SoA/Limiters.hpp"
#include <immintrin.h>
#include <cstring>
using namespace SoA;

Limiter::Limiter(const Domain& domain) : domain(domain) {
    Umin = (double*) _mm_malloc(sizeof(double)*(domain.Nc), ALIGN);
    Umax = (double*) _mm_malloc(sizeof(double)*(domain.Nc), ALIGN);
    dL   = (double*) _mm_malloc(sizeof(double)*(domain.Nf), ALIGN);
    dR   = (double*) _mm_malloc(sizeof(double)*(domain.Nf), ALIGN);
}

Limiter::~Limiter() {
    _mm_free(Umin);
    _mm_free(Umax);
    _mm_free(dL);
    _mm_free(dR);
}

void Limiter::limit_gradient(const int nvars, const double* const U[NVARS], const double* const G[NVARS], double* phi[NVARS]) {
    double eps = 1e-9;
    int Nf = domain.Nf;
    int Nc = domain.Nc;

    for (int var = 0; var < nvars; var++) {
        std::memcpy(Umin, U[var], Nc*sizeof(double));
        std::memcpy(Umax, U[var], Nc*sizeof(double));
        setarr(0, Nc, 1.0, phi[var]);

        for (int f = 0; f < Nf; f++) {
            int CL = domain.C[0][f];
            int CR = domain.C[1][f];
        
            Umin[CL] = fmin(Umin[CL], U[var][CR]);
            Umax[CL] = fmax(Umax[CL], U[var][CR]);
            Umin[CR] = fmin(Umin[CR], U[var][CL]);
            Umax[CR] = fmax(Umax[CR], U[var][CL]);
        }

        double rl[NDIMS], rr[NDIMS];
        double GL[NDIMS], GR[NDIMS];
        for (int f = 0; f < Nf; f++) {
            int CL = domain.C[0][f];
            int CR = domain.C[1][f];

            for (int dim = 0; dim < NDIMS; dim++) {
                rl[dim] = domain.rf[dim][0][f];
                rr[dim] = domain.rf[dim][1][f];
            }

            for (int dim = 0; dim < NDIMS; dim++) {
                GL[dim] = G[var][dim*Nc+CL];
                GR[dim] = G[var][dim*Nc+CR];
            }

            dL[f] = dot(GL, rl);
            dR[f] = dot(GR, rr);
        }

        for (int f = 0; f < Nf; f++) {
            int CL = domain.C[0][f];
            int CR = domain.C[1][f];

            double UL = U[var][CL];
            double UR = U[var][CR];

            double phi_fl = 1.0;
            double phi_fr = 1.0;

            double dLf = dL[f];
            double dRf = dR[f];

            phi_fl = (dLf >  eps) ? fmin(1, (Umax[CL]-UL)/dLf) : phi_fl;
            phi_fl = (dLf < -eps) ? fmin(1, (Umin[CL]-UL)/dLf) : phi_fl;

            phi_fr = (dRf >  eps) ? fmin(1, (Umax[CR]-UR)/dRf) : phi_fr;
            phi_fr = (dRf < -eps) ? fmin(1, (Umin[CR]-UR)/dRf) : phi_fr;

            phi[var][CL] = fmin(phi[var][CL], phi_fl);
            phi[var][CR] = fmin(phi[var][CR], phi_fr);
        }
    }
}
