#pragma once
#include "Global.hpp"
#include "AoSoA/Domain.hpp"
namespace AoSoA {

class Limiter {
    public:
    Limiter(const Domain& domain);
    ~Limiter();
    void limit_gradient(const int nvars, const double* const U[NVARS], const double* const G[NVARS], double* phi[NVARS]);

    private:
    const Domain& domain;
    double* Umin;
    double* Umax;
    double (*dU)[2][VECLEN];
};

}