#include "Base/Limiters.hpp"
#include <cstring>
using namespace Base;

Limiter::Limiter(const Domain& domain) : domain(domain) {
    Umin = new double[domain.Nc];
    Umax = new double[domain.Nc];
    dL   = new double[domain.Nf];
    dR   = new double[domain.Nf];
}

Limiter::~Limiter() {
    delete[] Umin;
    delete[] Umax;
    delete[] dL;
    delete[] dR;
}

void Limiter::limit_gradient(const int nvars, const double* const U[NVARS], const double* const G[NVARS], double* phi[NVARS]) {
    double eps = 1e-9;
    int Nf = domain.Nf;
    int Nc = domain.Nc;
    
    double GL[NDIMS], GR[NDIMS];

    for (int var = 0; var < nvars; var++) {
        std::memcpy(Umin, U[var], Nc*sizeof(double));
        std::memcpy(Umax, U[var], Nc*sizeof(double));
        setarr(0, Nc, 1.0, phi[var]);

        for (int f = 0; f < Nf; f++) {
            int CL = domain.C[f][0];
            int CR = domain.C[f][1];
        
            Umin[CL] = fmin(Umin[CL], U[var][CR]);
            Umax[CL] = fmax(Umax[CL], U[var][CR]);
            Umin[CR] = fmin(Umin[CR], U[var][CL]);
            Umax[CR] = fmax(Umax[CR], U[var][CL]);
        }

        for (int f = 0; f < Nf; f++) {
            int CL = domain.C[f][0];
            int CR = domain.C[f][1];

            for (int dim = 0; dim < NDIMS; dim++) {
                GL[dim] = G[var][dim*Nc+CL];
                GR[dim] = G[var][dim*Nc+CR];
            }

            dL[f] = dot(GL, domain.rf[f][0]);
            dR[f] = dot(GR, domain.rf[f][1]);
        }

        for (int f = 0; f < Nf; f++) {
            int CL = domain.C[f][0];
            int CR = domain.C[f][1];

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
