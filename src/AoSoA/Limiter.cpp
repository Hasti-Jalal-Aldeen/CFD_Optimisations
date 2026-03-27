#include "AoSoA/Limiters.hpp"
#include <cstring>
#include <immintrin.h>
using namespace AoSoA;

Limiter::Limiter(const Domain& domain) : domain(domain) {
    Umin = new double[domain.Nc];
    Umax = new double[domain.Nc];
    dU   = (double (*)[2][VECLEN]) _mm_malloc(sizeof(double)*domain.Nfb*2*VECLEN, ALIGN);
}

Limiter::~Limiter() {
    delete[] Umin;
    delete[] Umax;
    _mm_free(dU);
}

void Limiter::limit_gradient(const int nvars, const double* const U[NVARS], const double* const G[NVARS], double* phi[NVARS]) {
    double eps = 1e-9;
    int Nfb = domain.Nfb;
    int Nc  = domain.Nc;
    
    int *CL, *CR;

    double UL[VECLEN]        ALIGNED; 
    double UR[VECLEN]        ALIGNED; 
    double phi_fl[VECLEN]    ALIGNED;
    double phi_fr[VECLEN]    ALIGNED;
    double GL[NDIMS][VECLEN] ALIGNED; 
    double GR[NDIMS][VECLEN] ALIGNED;
    double rl[NDIMS][VECLEN] ALIGNED; 
    double rr[NDIMS][VECLEN] ALIGNED;

    for (int var = 0; var < nvars; var++) {
        std::memcpy(Umin, U[var], Nc*sizeof(double));
        std::memcpy(Umax, U[var], Nc*sizeof(double));
        setarr(0, Nc, 1.0, phi[var]);

        for (int fb = 0; fb < Nfb; fb++) {
            CL = domain.C[fb][0];
            CR = domain.C[fb][1];
        
            for (int v = 0; v < VECLEN; v++) {
                Umin[CL[v]] = fmin(Umin[CL[v]], U[var][CR[v]]);
                Umax[CL[v]] = fmax(Umax[CL[v]], U[var][CR[v]]);
                Umin[CR[v]] = fmin(Umin[CR[v]], U[var][CL[v]]);
                Umax[CR[v]] = fmax(Umax[CR[v]], U[var][CL[v]]);
            }
        }

        for (int fb = 0; fb < Nfb; fb++) {
            CL = domain.C[fb][0];
            CR = domain.C[fb][1];

            for (int dim = 0; dim < NDIMS; dim++) {
                for (int v = 0; v < VECLEN; v++) {
                    GL[dim][v] = G[var][dim*Nc+CL[v]];
                    GR[dim][v] = G[var][dim*Nc+CR[v]];
                    rl[dim][v] = domain.rf[fb][0][dim][v];
                    rr[dim][v] = domain.rf[fb][1][dim][v];
                }
            }

            dotv(GL, rl, dU[fb][0]);
            dotv(GR, rr, dU[fb][1]);
        }

        for (int fb = 0; fb < Nfb; fb++) {
            CL = domain.C[fb][0];
            CR = domain.C[fb][1];

            for (int v = 0; v < VECLEN; v++) {
                UL[v] = U[var][CL[v]];
                UR[v] = U[var][CR[v]];
            }

            for (int v = 0; v < VECLEN; v++) {
                phi_fl[v] = 1.0;
                phi_fr[v] = 1.0;
            }

            for (int v = 0; v < VECLEN; v++) {
                phi_fl[v] = (dU[fb][0][v] >  eps) ? fmin(1, (Umax[CL[v]]-UL[v])/dU[fb][0][v]) : phi_fl[v];
                phi_fl[v] = (dU[fb][0][v] < -eps) ? fmin(1, (Umin[CL[v]]-UL[v])/dU[fb][0][v]) : phi_fl[v];

                phi_fr[v] = (dU[fb][1][v] >  eps) ? fmin(1, (Umax[CR[v]]-UR[v])/dU[fb][1][v]) : phi_fr[v];
                phi_fr[v] = (dU[fb][1][v] < -eps) ? fmin(1, (Umin[CR[v]]-UR[v])/dU[fb][1][v]) : phi_fr[v];
            }

            for (int v = 0; v < VECLEN; v++) {
                phi[var][CL[v]] = fmin(phi[var][CL[v]], phi_fl[v]);
                phi[var][CR[v]] = fmin(phi[var][CR[v]], phi_fr[v]);
            }
        }
    }
}
