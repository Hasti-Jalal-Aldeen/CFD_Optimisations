#pragma once
#include "Global.hpp" 
#include "FAoS_CAoS/Domain.hpp"
namespace FAoS_CAoS {

class GradientScheme{
    public:
    GradientScheme(const Domain& domain);
    virtual void compute_gradient(int nvars, const double u[][NVARS], double grad_u[][NVARS][NDIMS]) const; 

    const int Nc, Nf;
    const Domain& domain;
};

class GaussGreen : public GradientScheme {
    public:
    using GradientScheme::GradientScheme;
    void compute_gradient(int nvars, const double u[][NVARS], double grad_u[][NVARS][NDIMS]) const override;
};

class NoGrads : public GradientScheme {
    public:
    using GradientScheme::GradientScheme;
};

}