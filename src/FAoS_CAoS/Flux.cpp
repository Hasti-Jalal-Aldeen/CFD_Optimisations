#include <cmath>
#include "FAoS_CAoS/Flux.hpp"
using namespace FAoS_CAoS;

Flux::Flux(int nvars, double gma, const Domain& domain) : nvars(nvars), gma(gma), domain(domain) {
    return;
}

void Flux::ComputeFlux(const double U[][NVARS], const double G[][NVARS][NDIMS], const double phi[][NVARS], double R[][NVARS]) const {
    const int Nf = domain.Nf;
    const int Nc = domain.Nc;

    double UfL[NVARS], UfR[NVARS], Ff[NVARS];

    for (int f = 0; f < Nf; f++) {
        int CL = domain.C[f][0];
        int CR = domain.C[f][1];
        double Af = domain.Af[f];

        MUSCL(U[CL], G[CL], phi[CL], domain.rf[f][0], UfL);
        MUSCL(U[CR], G[CR], phi[CR], domain.rf[f][1], UfR);
        Rusanov(UfL, UfR, domain.nf[f], Ff);

        for (int var = 0; var < nvars; var++) {
            R[CL][var] -= Af*Ff[var];
            R[CR][var] += Af*Ff[var];
        }
    }

    for (int var = 0; var < nvars; var++) {
        for (int i = 0; i < Nc; i++) { 
            R[i][var] *= domain.iv[i]; 
        }
    }
}

FORCE_INLINE void Flux::Rusanov(const double UL[NVARS], const double UR[NVARS], const double n[NDIMS], double F[NVARS]) const {
    double gmo = gma-1;

    double irl = 1.0/UL[0];
    double irr = 1.0/UR[0];

    double nvl = dot(n, &UL[1])*irl;
    double nvr = dot(n, &UR[1])*irr;

    double qql = 0.5*dot(&UL[1], &UL[1])*irl;
    double qqr = 0.5*dot(&UR[1], &UR[1])*irr;

    double pl  = gmo*(UL[NDIMS+1] - qql); 
    double pr  = gmo*(UR[NDIMS+1] - qqr); 

    double al  = sqrt(gma*pl*irl);
    double ar  = sqrt(gma*pr*irr);

    double S1  = fmax(fabs(nvl-al), fabs(nvr-ar));
    double S2  = fmax(fabs(nvl+al), fabs(nvr+ar));
    double S   = fmax(S1, S2);

    double FL[NVARS], FR[NVARS];
    EulerFlux(UL, n, FL);
    EulerFlux(UR, n, FR);

    for (int var = 0; var < NVARS; var++) { 
        F[var] = 0.5 * ( (FL[var] + FR[var]) - S*(UR[var] - UL[var]) );
    }
}

FORCE_INLINE void Flux::MUSCL(const double U[NVARS], const double G[NVARS][NDIMS], const double phi[NVARS], const double r[NDIMS], double Uf[NVARS]) const {
    for (int var = 0; var < NVARS; var++) { 
        Uf[var] = U[var] + phi[var]*dot(G[var], r);
    }
}

FORCE_INLINE void Flux::EulerFlux(const double U[NVARS], const double n[NDIMS], double F[NVARS]) const {
    double gmo = gma-1;
    double ir = 1.0/U[0];
    double P = gmo* (U[NDIMS+1] - 0.5*ir*dot(&U[1],&U[1]));
    #if NDIMS == 2
    F[0] = n[0]* U[1]              + n[1]* U[2]             ;
    F[1] = n[0]*(U[1]*U[1]*ir + P) + n[1]* U[2]*U[1]*ir     ;
    F[2] = n[0]* U[1]*U[2]*ir      + n[1]*(U[2]*U[2]*ir + P);
    F[3] = n[0]*(U[3]+P)*U[1]*ir   + n[1]*(U[3]+P)*U[2]*ir  ;
    #elif NDIMS == 3
    F[0] = n[0]* U[1]              + n[1]* U[2]              + n[2]* U[3]             ;
    F[1] = n[0]*(U[1]*U[1]*ir + P) + n[1]* U[2]*U[1]*ir      + n[2]* U[3]*U[1]*ir     ;
    F[2] = n[0]* U[1]*U[2]*ir      + n[1]*(U[2]*U[2]*ir + P) + n[2]* U[3]*U[2]*ir     ;
    F[3] = n[0]* U[1]*U[3]*ir      + n[1]* U[2]*U[3]*ir      + n[2]*(U[3]*U[3]*ir + P);
    F[4] = n[0]*(U[4]+P)*U[1]*ir   + n[1]*(U[4]+P)*U[2]*ir   + n[2]*(U[4]+P)*U[3]*ir  ;
    #endif
}



/*
void Flux::ComputeFlux(const double U[][NVARS], const double G[][NVARS][NDIMS], const double phi[][NVARS], double R[][NVARS]) const {
    const int Nf = domain.Nf;
    const int Nc = domain.Nc;

    int    *CL, *CR;
    double *Af;
    double nf   [NDIMS][VECLEN]        ALIGNED;
    double rL   [NDIMS][VECLEN]        ALIGNED;
    double rR   [NDIMS][VECLEN]        ALIGNED;

    double ULT  [NVARS]                ALIGNED;
    double URT  [NVARS]                ALIGNED;
    double phiLT[NVARS]                ALIGNED;
    double phiRT[NVARS]                ALIGNED;
    double GLT  [NVARS][NDIMS]         ALIGNED;
    double GRT  [NVARS][NDIMS]         ALIGNED;

    double UL   [NVARS][VECLEN]        ALIGNED;
    double UR   [NVARS][VECLEN]        ALIGNED;
    double phiL [NVARS][VECLEN]        ALIGNED;
    double phiR [NVARS][VECLEN]        ALIGNED;
    double GL   [NVARS][NDIMS][VECLEN] ALIGNED;
    double GR   [NVARS][NDIMS][VECLEN] ALIGNED;

    double UfL  [NVARS][VECLEN]        ALIGNED;
    double UfR  [NVARS][VECLEN]        ALIGNED;
    double Ff   [NVARS][VECLEN]        ALIGNED;

    for (int fb = 0; fb < Nf; fb++) { 
        CL = &domain.C[0][fb];
        CR = &domain.C[1][fb];

        for (int dim = 0; dim < NDIMS; dim++) {
            for (int v = 0; v < VECLEN; v++) {
                nf[dim][v] = domain.nf[dim][fb+v]   ;
                rL[dim][v] = domain.rf[dim][0][fb+v];
                rR[dim][v] = domain.rf[dim][1][fb+v];
            }
        }
        Af = &domain.Af[fb];

        for (int v = 0; v < VECLEN; v++) {
            #pragma omp simd
            for (int var = 0; var < NVARS; var++) {
                ULT[var]   = U[CL[v]][var];
                URT[var]   = U[CR[v]][var];
                phiLT[var] = phi[CL[v]][var];
                phiRT[var] = phi[CR[v]][var];
            }
            
            for (int var = 0; var < nvars; var++) {
                #pragma omp simd
                for (int dim = 0; dim < NDIMS; dim++) {
                    GLT[var][dim] = G[CL[v]][var][dim];
                    GRT[var][dim] = G[CR[v]][var][dim];
                }
            }

            #pragma omp simd
            for (int var = 0; var < NVARS; var++) {
                UL[var][v]   = ULT[var];
                UR[var][v]   = URT[var];
                phiL[var][v] = phiLT[var];
                phiR[var][v] = phiLT[var];
            }

            for (int var = 0; var < nvars; var++) {
                #pragma omp simd
                for (int dim = 0; dim < NDIMS; dim++) {
                    GL[var][dim][v] = GLT[var][dim];
                    GR[var][dim][v] = GRT[var][dim];
                }
            }

        }

        // MUSCL(UL, GL, phiL, rL, UfL);
        // MUSCL(UR, GR, phiR, rR, UfR);
        // Rusanov(UfL, UfR, nf, Ff);
        Rusanov(UL, UR, nf, Ff);

        for (int var = 0; var < nvars; var++) { 
            for (int v = 0; v < VECLEN; v++) {
                R[CL[v]][var] -= Af[v]*Ff[var][v];
                R[CR[v]][var] += Af[v]*Ff[var][v];
            }
        }
    }

    for (int i = 0; i < Nc; i++) { 
        for (int var = 0; var < NVARS; var++) { 
            R[i][var] *= domain.iv[i];
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

    for (int v = 0; v < VECLEN; v++) {
        irl[v] = 1.0/UL[0][v];
        irr[v] = 1.0/UR[0][v];
    }
    
    dotv(n, &UL[1], nvl);
    dotv(n, &UR[1], nvr);
    for (int v = 0; v < VECLEN; v++) {
        nvl[v] *= irl[v];
        nvr[v] *= irr[v];
    }

    dotv(&UL[1], &UL[1], qql);
    dotv(&UR[1], &UR[1], qqr);
    for (int v = 0; v < VECLEN; v++) {
        qql[v] *= 0.5*irl[v];
        qqr[v] *= 0.5*irr[v];
    }

    for (int v = 0; v < VECLEN; v++) {
        pl[v] = gmo*(UL[NDIMS+1][v] - qql[v]); 
        pr[v] = gmo*(UR[NDIMS+1][v] - qqr[v]); 
    }

    for (int v = 0; v < VECLEN; v++) {
        al[v] = sqrt(gma*pl[v]*irl[v]);
        ar[v] = sqrt(gma*pr[v]*irr[v]);
    }

    for (int v = 0; v < VECLEN; v++) {
        S1[v] = fmax(fabs(nvl[v]-al[v]), fabs(nvr[v]-ar[v]));
        S2[v] = fmax(fabs(nvl[v]+al[v]), fabs(nvr[v]+ar[v]));
    }
    for (int v = 0; v < VECLEN; v++) {
        S[v] = fmax(S1[v], S2[v]);
    }

    EulerFlux(UL, n, FL);
    EulerFlux(UR, n, FR);

    for (int var = 0; var < NVARS; var++) {
        for (int v = 0; v < VECLEN; v++) {
            F[var][v] = 0.5 * ( (FL[var][v] + FR[var][v]) - S[v]*(UR[var][v] - UL[var][v]) );
        }
    }
}

FORCE_INLINE void Flux::MUSCL(const double U[NVARS][VECLEN], const double G[NVARS][NDIMS][VECLEN], const double phi[NVARS][VECLEN], const double r[NDIMS][VECLEN], double Uf[NVARS][VECLEN]) const {
    double Gr[VECLEN] ALIGNED;
    for (int var = 0; var < NVARS; var++) {
        dotv(G[var], r, Gr);
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

    for (int v = 0; v < VECLEN; v++) {
        ir[v] = 1.0/U[0][v];
    }
    dotv(&U[1], &U[1], qq);
    for (int v = 0; v < VECLEN; v++) {
        P[v] = gmo* (U[NDIMS+1][v] - 0.5*ir[v]*qq[v]);
    }

    for (int v = 0; v < VECLEN; v++) {
        F[0][v] = n[0][v]* U[1][v]                       + n[1][v]* U[2][v]                      ;
        F[1][v] = n[0][v]*(U[1][v]*U[1][v]*ir[v] + P[v]) + n[1][v]* U[2][v]*U[1][v]*ir[v]        ;
        F[2][v] = n[0][v]* U[1][v]*U[2][v]*ir[v]         + n[1][v]*(U[2][v]*U[2][v]*ir[v] + P[v]);
        F[3][v] = n[0][v]*(U[3][v]+P[v])*U[1][v]*ir[v]   + n[1][v]*(U[3][v]+P[v])*U[2][v]*ir[v]  ;
    }
}
*/