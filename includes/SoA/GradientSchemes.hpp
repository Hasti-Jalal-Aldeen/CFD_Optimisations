#pragma once
#include "Global.hpp" 
#include "SoA/Domain.hpp"
namespace SoA {

class GradientScheme{
    public:
    GradientScheme(const Domain& domain);
    virtual void compute_gradient(int nvars, const double* const u[NVARS], double* grad_u[NVARS]) const; 

    const int Nc, Nf;
    const Domain& domain;
};

class GaussGreen : public GradientScheme {
    public:
    using GradientScheme::GradientScheme;
    void compute_gradient(int nvars, const double* const u[NVARS], double* grad_u[NVARS]) const override;
};

class NoGrads : public GradientScheme {
    public:
    using GradientScheme::GradientScheme;
};

}