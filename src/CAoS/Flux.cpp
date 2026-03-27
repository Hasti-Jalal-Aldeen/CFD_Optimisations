#include <cmath>
#include "CAoS/Flux.hpp"
using namespace CAoS;

Flux::Flux(int nvars, double gma, const Domain& domain) : nvars(nvars), gma(gma), domain(domain) {
    return;
}

void Flux::ComputeFlux(const double U[][NVARS], const double G[][NVARS][NDIMS], const double phi[][NVARS], double R[][NVARS]) const {
    const int Nf = domain.Nf;
    const int Nc = domain.Nc;

    double nf [NDIMS], rL [NDIMS], rR[NDIMS];
    double UfL[NVARS], UfR[NVARS], Ff[NVARS];

    for (int f = 0; f < Nf; f++) {
        int CL = domain.C[0][f];
        int CR = domain.C[1][f];

        for (int dim = 0; dim < NDIMS; dim++) {
            nf[dim] = domain.nf[dim][f];
            rL[dim] = domain.rf[dim][0][f];
            rR[dim] = domain.rf[dim][1][f];
        }
        double Af = domain.Af[f];

        MUSCL(U[CL], G[CL], phi[CL], rL, UfL);
        MUSCL(U[CR], G[CR], phi[CR], rR, UfR);
        Rusanov(UfL, UfR, nf, Ff);

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

