#pragma once
#include "Global.hpp"
#include "CAoS/Domain.hpp"
namespace CAoS {

class Limiter {
    public:
    Limiter(const Domain& domain);
    ~Limiter();
    void limit_gradient(const int nvars, const double U[][NVARS], const double G[][NVARS][NDIMS], double phi[][NVARS]);

    private:
    const Domain& domain;
    double (*Umin)[NVARS];
    double (*Umax)[NVARS];
    double (*dL)[NVARS];
    double (*dR)[NVARS]; 
};

}