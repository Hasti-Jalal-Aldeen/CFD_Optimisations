#pragma once
#include "Global.hpp"
#include "CAoS/Domain.hpp"
namespace CAoS {

enum BCtype {
    SYMMETRY ,
    WALL     ,
    SUPINLET ,
    SUBINLET ,
    SUPOUTLET,
    SUBOUTLET
};

class BoundaryCondition {
    public:
    BoundaryCondition(const Domain& domain);
    void init(int nvars, int F1, int F2, double gma, const double QBC[NVARS]);
    virtual void update_cells(double U[][NVARS]);
    void update_grads(double G[][NVARS][NDIMS]);
    void update_limiter(double phi[][NVARS]);
    
    int nvars, Nc;
    int F1, F2;
    double gma;
    double QBC[NVARS];
    double UBC[NVARS];
    const Domain& domain;
};

class Symmetry : public BoundaryCondition {
    public:
    using BoundaryCondition::BoundaryCondition;
    void update_cells(double U[][NVARS]) override;
};

class Wall : public BoundaryCondition {
    public:
    using BoundaryCondition::BoundaryCondition;
    void update_cells(double U[][NVARS]) override;
};

class SupInlet : public BoundaryCondition {
    public:
    using BoundaryCondition::BoundaryCondition;
    void update_cells(double U[][NVARS]) override;
};

class SupOutlet : public BoundaryCondition {
    public:
    using BoundaryCondition::BoundaryCondition;
    void update_cells(double U[][NVARS]) override;
};

class SubInlet : public BoundaryCondition {
    public:
    using BoundaryCondition::BoundaryCondition;
    void update_cells(double U[][NVARS]) override;
};

class SubOutlet : public BoundaryCondition {
    public:
    using BoundaryCondition::BoundaryCondition;
    void update_cells(double U[][NVARS]) override;
};

}