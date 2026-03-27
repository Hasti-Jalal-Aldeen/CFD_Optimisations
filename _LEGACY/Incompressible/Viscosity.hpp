#pragma once
#include "Global.hpp"
#include "Domain.hpp"

class ViscousModel {
    public:
    ViscousModel(double nu, const Domain& domain);
    virtual void discretise_viscosity(double* A) const;

    double nu;
    const Domain& domain;
};

class Inviscid : public ViscousModel {
    using ViscousModel::ViscousModel;
};

class DNS : public ViscousModel {
    public:
    using ViscousModel::ViscousModel;
    void discretise_viscosity(double* A) const override;
};
