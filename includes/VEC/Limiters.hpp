#pragma once
#include "Global.hpp"
#include "VEC/Domain.hpp"
namespace VEC {

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