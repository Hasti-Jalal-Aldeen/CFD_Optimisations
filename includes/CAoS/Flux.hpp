#pragma once
#include <immintrin.h>
#include "CAoS/Domain.hpp"
namespace CAoS {

class Flux {
    public:
    const int nvars;
    const double gma;
    Flux(int nvars, double gma, const Domain& domain);

    void ComputeFlux(const double U[][NVARS], const double G[][NVARS][NDIMS], const double phi[][NVARS], double R[][NVARS]) const;

    private:
    const Domain& domain;

    FORCE_INLINE void EulerFlux(const double U[NVARS], const double n[NDIMS], double F[NVARS]) const;
    FORCE_INLINE void Rusanov(const double UL[NVARS], const double UR[NVARS], const double n[NDIMS], double F[NVARS]) const;
    FORCE_INLINE void MUSCL(const double U[NVARS], const double G[NVARS][NDIMS], const double phi[NVARS], const double r[NDIMS], double Uf[NVARS]) const;
};

}
