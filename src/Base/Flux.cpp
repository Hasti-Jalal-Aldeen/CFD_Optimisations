#include <cmath>
#include "Base/Flux.hpp"
using namespace Base;

Flux::Flux(int nvars, double gma, const Domain& domain) : nvars(nvars), gma(gma), domain(domain) {
    return;
}

void Flux::ComputeFlux(const double* const U[NVARS], const double* const G[NVARS], const double* const phi[NVARS], double* R[NVARS]) const {
    const int Nf = domain.Nf;
    const int Nc = domain.Nc;

    double UL  [NVARS]       , UR  [NVARS]       ;
    double phiL[NVARS]       , phiR[NVARS]       ;
    double GL  [NVARS][NDIMS], GR  [NVARS][NDIMS];

    double UfL[NVARS], UfR[NVARS], Ff[NVARS];

    for (int f = 0; f < Nf; f++) {
        int CL = domain.C[f][0];
        int CR = domain.C[f][1];
        double Af = domain.Af[f];

        for (int var = 0; var < nvars; var++) { 
            UL  [var] = U  [var][CL];
            UR  [var] = U  [var][CR];
            phiL[var] = phi[var][CL];
            phiR[var] = phi[var][CR];
            for (int dim = 0; dim < NDIMS; dim++) {
                GL[var][dim] = G[var][dim*Nc+CL];
                GR[var][dim] = G[var][dim*Nc+CR];
            }
        }

        MUSCL(UL, GL, phiL, domain.rf[f][0], UfL);
        MUSCL(UR, GR, phiR, domain.rf[f][1], UfR);
        Rusanov(UfL, UfR, domain.nf[f], Ff);

        for (int var = 0; var < nvars; var++) {
            R[var][CL] -= Af*Ff[var];
            R[var][CR] += Af*Ff[var];
        }
    }

    for (int var = 0; var < nvars; var++) {
        for (int i = 0; i < Nc; i++) { 
            R[var][i] *= domain.iv[i]; 
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
    // double un = dot(n, &U[1])*ir;
    // for (int var = NDIMS+2; var < nvars; var++) {
    //     F[var] = un*U[var];
    // }
}

/*
void ViscFlux(const double U[4], const double G[4][2], const double n[2], double F[4], double gma, double mu, double Pr) {
    double r  = U[0]; 
    double ru = U[1];
    double rv = U[2];
    double rE = U[3];

    double rcprho = 1.0/r;
    double u = ru*rcprho;
    double v = rv*rcprho;

    double r_x = G[0][0];
    double r_y = G[0][1];

    // Vel grad (rho*grad[u,v])
    double u_x = G[1][0] - u*r_x;
    double u_y = G[1][1] - u*r_y;
    double v_x = G[2][0] - v*r_x;
    double v_y = G[2][1] - v*r_y;

    double rE_x = G[3][0];
    double rE_y = G[3][1];

    // Temp derivatives Cv*dT/d[x,y]
    double T_x = rcprho*(rE_x - (rcprho*r_x*rE + u*u_x + v*v_x));
    double T_y = rcprho*(rE_y - (rcprho*r_y*rE + u*u_y + v*v_y));

    double t_xx = 2*mu*rcprho*(u_x - (u_x + v_y)/3.0);
    double t_yy = 2*mu*rcprho*(v_y - (u_x + v_y)/3.0);
    double t_xy = mu*rcprho*(v_x + u_y);

    double FE_x = u*t_xx + v*t_xy - mu*gma*T_x/Pr;
    double FE_y = u*t_xy + v*t_yy - mu*gma*T_y/Pr;

    // F[0] -= 0;
    F[1] -= n[0]*t_xx + n[1]*t_xy;
    F[2] -= n[0]*t_xy + n[1]*t_yy;
    F[3] -= n[0]*FE_x + n[1]*FE_y;
}
void ViscFlux(const double U[5], const double G[5][3], const double n[3], double F[5], double gma, double mu, double Pr) {
    double r = U[0]; 
    double ru = U[1];
    double rv = U[2];
    double rw = U[3];
    double rE = U[4];

    double rcprho = 1.0/r;
    double u = ru*rcprho;
    double v = rv*rcprho;
    double w = rw*rcprho;

    double r_x = G[0][0];
    double r_y = G[0][1];
    double r_z = G[0][2];

    // Vel grad (rho*grad[u,v])
    double u_x = G[1][0] - u*r_x;
    double u_y = G[1][1] - u*r_y;
    double u_z = G[1][2] - u*r_z;

    double v_x = G[2][0] - v*r_x;
    double v_y = G[2][1] - v*r_y;
    double v_z = G[2][2] - v*r_z;

    double w_x = G[3][0] - w*r_x;
    double w_y = G[3][1] - w*r_y;
    double w_z = G[3][2] - w*r_z;

    double rE_x = G[4][0];
    double rE_y = G[4][1];
    double rE_z = G[4][2];

    // Temp derivatives Cv*dT/d[x,y]
    double T_x = rcprho*(rE_x - (rcprho*r_x*rE + u*u_x + v*v_x + w*w_x));
    double T_y = rcprho*(rE_y - (rcprho*r_y*rE + u*u_y + v*v_y + w*w_y));
    double T_z = rcprho*(rE_z - (rcprho*r_z*rE + u*u_z + v*v_z + w*w_z));

    double t_xx = 2*mu*rcprho*(u_x - (u_x + v_y + w_z)/3.0);
    double t_yy = 2*mu*rcprho*(v_y - (u_x + v_y + w_z)/3.0);
    double t_zz = 2*mu*rcprho*(w_z - (u_x + v_y + w_z)/3.0);
    double t_xy = mu*rcprho*(v_x + u_y);
    double t_xz = mu*rcprho*(w_x + u_z);
    double t_yz = mu*rcprho*(w_y + v_z);

    double FE_x = u*t_xx + v*t_xy + w*t_xz - mu*gma*T_x/Pr;
    double FE_y = u*t_xy + v*t_yy + w*t_yz - mu*gma*T_y/Pr;
    double FE_z = u*t_xz + v*t_yz + w*t_zz - mu*gma*T_z/Pr;

    // F[0] -= 0;
    F[1] -= n[0]*t_xx + n[1]*t_xy + n[2]*t_xz;
    F[2] -= n[0]*t_xy + n[1]*t_yy + n[2]*t_yz;
    F[3] -= n[0]*t_xz + n[1]*t_yz + n[2]*t_zz;
    F[4] -= n[0]*FE_x + n[1]*FE_y + n[2]*FE_z;
}
*/