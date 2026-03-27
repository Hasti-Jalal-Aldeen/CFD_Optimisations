#pragma once
#include "Global.hpp"
#include "SoA/Domain.hpp"
namespace SoA {

class Limiter {
    public:
    Limiter(const Domain& domain);
    ~Limiter();
    void limit_gradient(const int nvars, const double* const U[NVARS], const double* const G[NVARS], double* phi[NVARS]);

    private:
    const Domain& domain;
    double* Umin;
    double* Umax;
    double* dL;
    double* dR; 
};

}