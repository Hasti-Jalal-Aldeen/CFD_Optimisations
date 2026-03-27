#pragma once
#include "Global.hpp"
#include "Domain.hpp"

class ConvectionScheme {
    public:
    ConvectionScheme();
    virtual void discretise_momentum(double* A, double* b[NDIMS], double* F, double* G[NDIMS], const Domain& domain) const;
};

class UD : public ConvectionScheme {
    public:
    using ConvectionScheme::ConvectionScheme;
    void discretise_momentum(double* A, double* b[NDIMS], double* F, double* G[NDIMS], const Domain& domain) const override;
};

class LUD : public ConvectionScheme {
    public:
    using ConvectionScheme::ConvectionScheme;
    void discretise_momentum(double* A, double* b[NDIMS], double* F, double* G[NDIMS], const Domain& domain) const override;
};
