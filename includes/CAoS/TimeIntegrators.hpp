#pragma once
#include "Global.hpp"
#include "CAoS/Flux.hpp"
#include "CAoS/GradientSchemes.hpp"
#include "CAoS/Domain.hpp"
#include "CAoS/Limiters.hpp"
#include "Writer.hpp"
#include "CAoS/BoundaryConditions.hpp"
namespace CAoS {

class Euler {
    public:
    Euler(int nvars, Domain& domain, const Flux& flux, const GradientScheme& grads, Limiter& limiter);
    ~Euler();
    void gather_BCs(const BCtype* BCmap, const double QBC[][NVARS]);
    void initialise(const double QIC[NVARS]);

    void run(double T, double dt, double dt_out);
    void compute_prims();

    private:
    const int Nc, Nf, nvars;
    Domain& domain;
    const Flux& flux;
    const GradientScheme& grads;
    Limiter& limiter;

    void add_residuals(double dt);

    // BCs
    const int NBC;
    BoundaryCondition** BCs;
    void update_BC_cells();
    void update_BC_grads();
    void update_BC_limiter();

    // Cells
    double* Q[NVARS]; // Primitive Variables (r, u, v, w, p)
    double (*U)[NVARS]; // Conservative Variables
    double (*R)[NVARS]; // rhs
    double (*phi)[NVARS]; // slope limiter
    double (*G)[NVARS][NDIMS]; // Conservative Gradients SoA

    Writer* writer;
};

}
