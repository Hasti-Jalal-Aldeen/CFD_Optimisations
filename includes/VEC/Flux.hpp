#pragma once
#include "VEC/Domain.hpp"
namespace VEC {

class Flux {
    public:
    const int nvars;
    const double gma;
    Flux(int nvars, double gma, const Domain& domain);

    void ComputeFlux(const double* const U[NVARS], const double* const G[NVARS], const double* const phi[NVARS], double* R[NVARS]) const;

    private:
    const Domain& domain;

    FORCE_INLINE void EulerFlux(const double U[NVARS][VECLEN], const double n[NDIMS][VECLEN], double F[NVARS][VECLEN]) const;
    FORCE_INLINE void Rusanov(const double UL[NVARS][VECLEN], const double UR[NVARS][VECLEN], const double n[NDIMS][VECLEN], double F[NVARS][VECLEN]) const;
    FORCE_INLINE void MUSCL(const double U[NVARS][VECLEN], const double G[NVARS][NDIMS][VECLEN], const double phi[NVARS][VECLEN], const double r[NDIMS][VECLEN], double Uf[NVARS][VECLEN]) const;
};

}
