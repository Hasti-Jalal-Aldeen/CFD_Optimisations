#include <cmath>
#include "VEC/Flux.hpp"
using namespace VEC;

Flux::Flux(int nvars, double gma, const Domain& domain) : nvars(nvars), gma(gma), domain(domain) {
    return;
}

void Flux::ComputeFlux(const double* const U[NVARS], const double* const G[NVARS], const double* const phi[NVARS], double* R[NVARS]) const {
    const int Nf = domain.Nf;
    const int Nc = domain.Nc;

    int    *CL, *CR;
    double *Af;
    double nf  [NDIMS][VECLEN]        ALIGNED;
    double rL  [NDIMS][VECLEN]        ALIGNED;
    double rR  [NDIMS][VECLEN]        ALIGNED;
    
    double UL  [NVARS][VECLEN]        ALIGNED;
    double UR  [NVARS][VECLEN]        ALIGNED;
    double phiL[NVARS][VECLEN]        ALIGNED;
    double phiR[NVARS][VECLEN]        ALIGNED;
    double GL  [NVARS][NDIMS][VECLEN] ALIGNED;
    double GR  [NVARS][NDIMS][VECLEN] ALIGNED;

    double UfL [NVARS][VECLEN]        ALIGNED;
    double UfR [NVARS][VECLEN]        ALIGNED;
    double Ff  [NVARS][VECLEN]        ALIGNED;

    for (int fb = 0; fb < Nf; fb+=VECLEN) { 
        CL = &domain.C[0][fb];
        CR = &domain.C[1][fb];

        for (int dim = 0; dim < NDIMS; dim++) {
            #pragma omp simd simdlen(VECLEN) safelen(VECLEN)
            for (int v = 0; v < VECLEN; v++) {
                nf[dim][v] = domain.nf[dim][fb+v]   ;
                rL[dim][v] = domain.rf[dim][0][fb+v];
                rR[dim][v] = domain.rf[dim][1][fb+v];
            }
        }
        Af = &domain.Af[fb];
  
        for (int var = 0; var < nvars; var++) {   
            #pragma omp simd simdlen(VECLEN) safelen(VECLEN)
            for (int v = 0; v < VECLEN; v++) {
                UL  [var][v] = U  [var][CL[v]];
                UR  [var][v] = U  [var][CR[v]];
                phiL[var][v] = phi[var][CL[v]];
                phiR[var][v] = phi[var][CR[v]];
            }
            for (int dim = 0; dim < NDIMS; dim++) {
                #pragma omp simd simdlen(VECLEN) safelen(VECLEN)
                for (int v = 0; v < VECLEN; v++) {
                    GL[var][dim][v] = G[var][dim*Nc+CL[v]];
                    GR[var][dim][v] = G[var][dim*Nc+CR[v]];
                }
            }
        }

        MUSCL(UL, GL, phiL, rL, UfL);
        MUSCL(UR, GR, phiR, rR, UfR);
        Rusanov(UfL, UfR, nf, Ff);

        for (int var = 0; var < nvars; var++) { 
            #pragma omp simd simdlen(VECLEN) safelen(VECLEN)
            for (int v = 0; v < VECLEN; v++) {
                R[var][CL[v]] -= Af[v]*Ff[var][v]; 
                R[var][CR[v]] += Af[v]*Ff[var][v];
            }
        }
    }

    for (int var = 0; var < nvars; var++) { 
        #pragma omp simd
        for (int i = 0; i < Nc; i++) { 
            R[var][i] *= domain.iv[i];
        }
    }
}

FORCE_INLINE void Flux::Rusanov(const double UL[NVARS][VECLEN], const double UR[NVARS][VECLEN], const double n[NDIMS][VECLEN], double F[NVARS][VECLEN]) const {
    double gmo = gma-1;

    double irl[VECLEN] ALIGNED;
    double irr[VECLEN] ALIGNED;

    double nvl[VECLEN] ALIGNED;
    double nvr[VECLEN] ALIGNED;

    double qql[VECLEN] ALIGNED;
    double qqr[VECLEN] ALIGNED;

    double pl [VECLEN] ALIGNED;
    double pr [VECLEN] ALIGNED;

    double al [VECLEN] ALIGNED;
    double ar [VECLEN] ALIGNED;

    double S1 [VECLEN] ALIGNED;
    double S2 [VECLEN] ALIGNED;
    double S  [VECLEN] ALIGNED;

    double FL [NVARS][VECLEN] ALIGNED; 
    double FR [NVARS][VECLEN] ALIGNED;
    
    #pragma omp simd simdlen(VECLEN) safelen(VECLEN)
    for (int v = 0; v < VECLEN; v++) {
        irl[v] = 1.0/UL[0][v];
        irr[v] = 1.0/UR[0][v];
    }
    
    dotv(n, &UL[1], nvl);
    dotv(n, &UR[1], nvr);
    #pragma omp simd simdlen(VECLEN) safelen(VECLEN)
    for (int v = 0; v < VECLEN; v++) {
        nvl[v] *= irl[v];
        nvr[v] *= irr[v];
    }

    dotv(&UL[1], &UL[1], qql);
    dotv(&UR[1], &UR[1], qqr);
    #pragma omp simd simdlen(VECLEN) safelen(VECLEN)
    for (int v = 0; v < VECLEN; v++) {
        qql[v] *= 0.5*irl[v];
        qqr[v] *= 0.5*irr[v];
    }

    #pragma omp simd simdlen(VECLEN) safelen(VECLEN)
    for (int v = 0; v < VECLEN; v++) {
        pl[v] = gmo*(UL[NDIMS+1][v] - qql[v]); 
        pr[v] = gmo*(UR[NDIMS+1][v] - qqr[v]); 
    }

    #pragma omp simd simdlen(VECLEN) safelen(VECLEN)
    for (int v = 0; v < VECLEN; v++) {
        al[v] = sqrt(gma*pl[v]*irl[v]);
        ar[v] = sqrt(gma*pr[v]*irr[v]);
    }

    #pragma omp simd simdlen(VECLEN) safelen(VECLEN)
    for (int v = 0; v < VECLEN; v++) {
        S1[v] = fmax(fabs(nvl[v]-al[v]), fabs(nvr[v]-ar[v]));
        S2[v] = fmax(fabs(nvl[v]+al[v]), fabs(nvr[v]+ar[v]));
    }
    #pragma omp simd simdlen(VECLEN) safelen(VECLEN)
    for (int v = 0; v < VECLEN; v++) {
        S[v] = fmax(S1[v], S2[v]);
    }

    EulerFlux(UL, n, FL);
    EulerFlux(UR, n, FR);

    for (int var = 0; var < NVARS; var++) { 
        #pragma omp simd simdlen(VECLEN) safelen(VECLEN)
        for (int v = 0; v < VECLEN; v++) {
            F[var][v] = 0.5 * ( (FL[var][v] + FR[var][v]) - S[v]*(UR[var][v] - UL[var][v]) );
        }
    }
}

FORCE_INLINE void Flux::MUSCL(const double U[NVARS][VECLEN], const double G[NVARS][NDIMS][VECLEN], const double phi[NVARS][VECLEN], const double r[NDIMS][VECLEN], double Uf[NVARS][VECLEN]) const {
    double Gr[VECLEN] ALIGNED;
    for (int var = 0; var < NVARS; var++) {  
        dotv(G[var], r, Gr);
        #pragma omp simd simdlen(VECLEN) safelen(VECLEN)
        for (int v = 0; v < VECLEN; v++) {
            Uf[var][v] = U[var][v] + phi[var][v]*Gr[v];
        }
    }
}

FORCE_INLINE void Flux::EulerFlux(const double U[NVARS][VECLEN], const double n[NDIMS][VECLEN], double F[NVARS][VECLEN]) const {
    double gmo = gma-1;
    double ir[VECLEN] ALIGNED;
    double qq[VECLEN] ALIGNED;
    double P [VECLEN] ALIGNED;

    #pragma omp simd simdlen(VECLEN) safelen(VECLEN)
    for (int v = 0; v < VECLEN; v++) {
        ir[v] = 1.0/U[0][v];
    }
    dotv(&U[1], &U[1], qq);
    #pragma omp simd simdlen(VECLEN) safelen(VECLEN)
    for (int v = 0; v < VECLEN; v++) {
        P[v] = gmo* (U[NDIMS+1][v] - 0.5*ir[v]*qq[v]);
    }

    #pragma omp simd simdlen(VECLEN) safelen(VECLEN)
    for (int v = 0; v < VECLEN; v++) {
        F[0][v] = n[0][v]* U[1][v]                       + n[1][v]* U[2][v]                      ;
        F[1][v] = n[0][v]*(U[1][v]*U[1][v]*ir[v] + P[v]) + n[1][v]* U[2][v]*U[1][v]*ir[v]        ;
        F[2][v] = n[0][v]* U[1][v]*U[2][v]*ir[v]         + n[1][v]*(U[2][v]*U[2][v]*ir[v] + P[v]);
        F[3][v] = n[0][v]*(U[3][v]+P[v])*U[1][v]*ir[v]   + n[1][v]*(U[3][v]+P[v])*U[2][v]*ir[v]  ;
    }
}
