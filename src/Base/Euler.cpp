#include "Base/TimeIntegrators.hpp"
using namespace Base;

Euler::Euler(int nvars, Domain& domain, const Flux& flux, const GradientScheme& grads, Limiter& limiter)
: nvars(nvars), Nc(domain.Nc), Nf(domain.Nf), NBC(domain.NBC), domain(domain), flux(flux), grads(grads), limiter(limiter)
{
    BCs = new BoundaryCondition*[NBC];

    for (int var = 0; var < nvars; var++) {
        Q[var]    = new double[Nc]();
        U[var]    = new double[Nc]();
        R[var]    = new double[Nc]();
        phi[var]  = new double[Nc]();
        G[var]    = new double[NDIMS*Nc]();
    }

    #if   NDIMS == 2
    std::string varnames[] = {"Density", "Velocity_X", "Velocity_Y", "Pressure"};
    #elif NDIMS == 3
    std::string varnames[] = {"Density", "Velocity_X", "Velocity_Y", "Velocity_Z", "Pressure"};
    #endif
    writer = new Writer(nvars, "../sols/Base", varnames, Nc, Q);

    domain.zero_bc_iv();
}

Euler::~Euler() {
    for (int bc = 0; bc < NBC; bc++) delete BCs[bc];
    delete[] BCs;

    for (int var = 0; var < nvars; var++) {
        delete[] Q[var];
        delete[] U[var];
        delete[] R[var];
        delete[] phi[var];
        delete[] G[var];
    }

    delete writer;
}

void Euler::gather_BCs(const BCtype* BCmap, const double QBC[][NVARS]) {
    int fcolBC = domain.Ncol-NBC;
    for (int fcol = fcolBC; fcol < domain.Ncol; fcol++) {
        int F1 = domain.mesh.coloff[fcol];
        int F2 = domain.mesh.coloff[fcol+1];
        int bc = fcol-fcolBC;

        BCtype BCtype = BCmap[ domain.mesh.coltag[fcol] ];
        switch (BCtype) {
            case SYMMETRY : BCs[bc] = new Symmetry(domain) ; break;
            case WALL     : BCs[bc] = new Wall(domain)     ; break;
            case SUPINLET : BCs[bc] = new SupInlet(domain) ; break;
            case SUBINLET : BCs[bc] = new SubInlet(domain) ; break;
            case SUPOUTLET: BCs[bc] = new SupOutlet(domain); break;
            case SUBOUTLET: BCs[bc] = new SubOutlet(domain); break;
        }
        BCs[bc]->init(nvars, F1, F2, flux.gma, QBC[bc]);
    }
}

void Euler::initialise(const double QIC[NVARS]) {
    double UIC[NVARS];
    UIC[0] = QIC[0];
    for (int var = 1; var < nvars; var++) {
        UIC[var] = QIC[0]*QIC[var];
    }
    UIC[NDIMS+1] = QIC[NDIMS+1]/(flux.gma-1) + 0.5*QIC[0]*dot(&QIC[1], &QIC[1]);

    for (int var = 0; var < nvars; var++) {
        for (int i = 0; i < Nc; i++) {
            U[var][i] = UIC[var];
        }
    }
    
    update_BC_cells();
}

void Euler::run(double T, double dt, double dt_out) {
    double t = 0.0;
    double t_out = dt_out;

    compute_prims();
    writer->write_soln();
    while (t < T) {
        grads.compute_gradient(nvars, U, G);
        limiter.limit_gradient(nvars, U, G, phi); 
        update_BC_grads();
        update_BC_limiter();

        flux.ComputeFlux(U, G, phi, R);
        add_residuals(dt);
        for (int var = 0; var < nvars; var++) {
            setarr(0, Nc, 0.0, R[var]); // clear residuals
        }
        update_BC_cells();
        t += dt;

        if (t >= t_out) {
            compute_prims();
            writer->write_soln();
            t_out += dt_out;
        }
    }
}

void Euler::add_residuals(double dt) {
    for (int var = 0; var < nvars; var++) {
        for (int i = 0; i < Nc; i++) {
            U[var][i] += dt*R[var][i];
        }
    }
}

void Euler::compute_prims() {
    for (int i = 0; i < Nc; i++) {
        Q[0][i] = U[0][i];
    }

    for (int var = 1; var < nvars; var++) {
        if (var == NDIMS+1) continue;
        for (int i = 0; i < Nc; i++) {
            Q[var][i] = U[var][i]/U[0][i];
        }
    }

    double gmo = flux.gma-1;
    for (int i = 0; i < Nc; i++) {
        double u[NDIMS];
        for (int dim = 0; dim < NDIMS; dim++) {
            u[dim] = Q[dim+1][i];
        }

        Q[NDIMS+1][i] = gmo * (U[NDIMS+1][i] - 0.5*Q[0][i]*dot(u,u));
    }
}

void Euler::update_BC_cells() {
    for (int bc = 0; bc < NBC; bc++) {
        BCs[bc]->update_cells(U);
    }
}

void Euler::update_BC_grads() {
    for (int bc = 0; bc < NBC; bc++) {
        BCs[bc]->update_grads(G);
    }
}

void Euler::update_BC_limiter() {
    for (int bc = 0; bc < NBC; bc++) {
        BCs[bc]->update_limiter(phi);
    }
}
