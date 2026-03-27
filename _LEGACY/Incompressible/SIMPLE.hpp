#pragma once
#include "Global.hpp"
#include "Domain.hpp"
#include "GradientSchemes.hpp"
#include "IterativeSolvers.hpp"
#include "Writer.hpp"
#include "Incompressible/BoundaryConditions.hpp"
#include "Incompressible/ConvectionSchemes.hpp"
#include "Incompressible/Viscosity.hpp"

class SIMPLE {
    public:
    SIMPLE(const Domain& domain, const GradientScheme& Ugrads, const GradientScheme& Pgrads, 
        const IterativeSolver& mom_solver, const IterativeSolver& pres_solver, const ConvectionScheme& convScheme, const ViscousModel& viscScheme);
    ~SIMPLE();
    void gather_BCs(const int* BCmap, const double UBC[][NDIMS], const double PBC[]);
    void initialise(double P_ref, const double UIC[NDIMS], double PIC);
    void run_solver(int niters);

    private:
    int iter;

    double* writer_arrs[NDIMS+1];
    Writer* writer;

    const int Nc, Nf;
    const Domain& domain;

    const ConvectionScheme& convScheme;
    const ViscousModel& viscScheme;
    const IterativeSolver& mom_solver;
    const IterativeSolver& pres_solver;
    void compute_F();
    void underrelax();
    void scale_mom_sys_by_vol();
    void solve_mom_eqn();
    void solve_pres_eqn();

    // BCs
    const int NBC;
    BoundaryCondition** BCs;
    void update_BC_cells();
    void update_BC_mom_matrix();
    void update_BC_pressure_matrix();
    
    // Grads
    const GradientScheme& Pgrads;
    const GradientScheme& Ugrads;

    // cells
    double P_ref;
    double* U[NDIMS];
    double* G[NDIMS];   // G[var][dim*Nc+i]
    double* P;          // Gauge pressure
    double* gradP;

    // momentum equation
    double* A_mom;
    double* D_mom; // diagonal of momentum matrix
    double* b_mom[NDIMS];
    double* rhs_mom[NDIMS];
    double* F;

    // pressure equation
    double* A_pres;
    double* H_a[NDIMS];
    double* rhs_pres;
    double* Pp;
};
