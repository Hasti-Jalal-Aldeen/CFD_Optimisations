#include "CAoS/Limiters.hpp"
#include <cstring>
using namespace CAoS;

Limiter::Limiter(const Domain& domain) : domain(domain) {
    Umin = new double[domain.Nc][NVARS];
    Umax = new double[domain.Nc][NVARS];
    dL   = new double[domain.Nf][NVARS];
    dR   = new double[domain.Nf][NVARS];
}

Limiter::~Limiter() {
    delete[] Umin;
    delete[] Umax;
    delete[] dL;
    delete[] dR;
}

void Limiter::limit_gradient(const int nvars, const double U[][NVARS], const double G[][NVARS][NDIMS], double phi[][NVARS]) {
    double eps = 1e-9;
    int Nf = domain.Nf;
    int Nc = domain.Nc;
    
    std::memcpy(Umin[0], U[0], NVARS*Nc*sizeof(double));
    std::memcpy(Umax[0], U[0], NVARS*Nc*sizeof(double));
    setarr(0, NVARS*Nc, 1.0, phi[0]);

    for (int f = 0; f < Nf; f++) {
        int CL = domain.C[0][f];
        int CR = domain.C[1][f];
        
        for (int var = 0; var < NVARS; var++) {
            Umin[CL][var] = fmin(Umin[CL][var], U[CR][var]);
            Umax[CL][var] = fmax(Umax[CL][var], U[CR][var]);
            Umin[CR][var] = fmin(Umin[CR][var], U[CL][var]);
            Umax[CR][var] = fmax(Umax[CR][var], U[CL][var]);
        }
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

        for (int var = 0; var < NVARS; var++) {
            for (int dim = 0; dim < NDIMS; dim++) {
                GL[dim] = G[CL][var][dim];
                GR[dim] = G[CR][var][dim];
            }
    
            dL[f][var] = dot(GL, rl);
            dR[f][var] = dot(GR, rr);
        }
    }
    
    double phi_fl[NVARS], phi_fr[NVARS];
    for (int f = 0; f < Nf; f++) {
        int CL = domain.C[0][f];
        int CR = domain.C[1][f];

        const double* UL = U[CL];
        const double* UR = U[CR];
        const double* dLf = dL[f];
        const double* dRf = dR[f];

        for (int var = 0; var < NVARS; var++) {
            phi_fl[var] = 1.0;
            phi_fr[var] = 1.0;
        }

        for (int var = 0; var < NVARS; var++) {
            phi_fl[var] = (dLf[var] >  eps) ? fmin(1, (Umax[CL][var]-UL[var])/dLf[var]) : phi_fl[var];
            phi_fl[var] = (dLf[var] < -eps) ? fmin(1, (Umin[CL][var]-UL[var])/dLf[var]) : phi_fl[var];
    
            phi_fr[var] = (dRf[var] >  eps) ? fmin(1, (Umax[CR][var]-UR[var])/dRf[var]) : phi_fr[var];
            phi_fr[var] = (dRf[var] < -eps) ? fmin(1, (Umin[CR][var]-UR[var])/dRf[var]) : phi_fr[var];
        }
        for (int var = 0; var < NVARS; var++) {
            phi[CL][var] = fmin(phi[CL][var], phi_fl[var]);
            phi[CR][var] = fmin(phi[CR][var], phi_fr[var]);
        }
    }
}
