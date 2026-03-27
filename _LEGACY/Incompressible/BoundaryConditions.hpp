#pragma once
#include "Global.hpp"
#include "Domain.hpp"

// Neumann_0 is a bc where there is 0 gradient at the BC
void modify_matrix_Dirichlet(int F1, int F2, int N, double* A, const Domain& domain);
void modify_matrix_Neumann  (int F1, int F2, int N, double* A, const Domain& domain);
void modify_rhs_Dirichlet   (int F1, int F2, double* rhs, double UBC  , const Domain& domain);
void modify_rhs_Neumann_0   (int F1, int F2, double* rhs,               const Domain& domain);
void modify_rhs_Neumann     (int F1, int F2, double* rhs, double du_dn, const Domain& domain);

class BoundaryCondition {
    public:
    BoundaryCondition(const Domain& domain);
    void init(int nvars, int F1, int F2, const double UBC[], double P_BC, double* U[NDIMS], double* P);
    virtual void modify_mom_matrix(int N, double* A, double* rhs[NDIMS]);
    virtual void modify_pressure_matrix(int N, double* A, double* rhs);
    virtual void update_cells();
    
    int nvars;
    int F1, F2;
    double UBC[NDIMS];
    double P_BC;

    const Domain& domain;

    double* U[NDIMS];
    double* P;
};

class Symmetry : public BoundaryCondition {
    public:
    using BoundaryCondition::BoundaryCondition;
    void modify_mom_matrix(int N, double* A, double* rhs[NDIMS]) override;
    void modify_pressure_matrix(int N, double* A, double* rhs) override;
    void update_cells() override;
};

class Wall : public BoundaryCondition {
    public:
    using BoundaryCondition::BoundaryCondition;
    void modify_mom_matrix(int N, double* A, double* rhs[NDIMS]) override;
    void modify_pressure_matrix(int N, double* A, double* rhs) override;
    void update_cells() override;
};

class Inlet : public BoundaryCondition {
    public:
    using BoundaryCondition::BoundaryCondition;
    void modify_mom_matrix(int N, double* A, double* rhs[NDIMS]) override;
    void modify_pressure_matrix(int N, double* A, double* rhs) override;
    void update_cells() override;
};

class Outlet : public BoundaryCondition {
    public:
    using BoundaryCondition::BoundaryCondition;
    void modify_mom_matrix(int N, double* A, double* rhs[NDIMS]) override;
    void modify_pressure_matrix(int N, double* A, double* rhs) override;
    void update_cells() override;
};