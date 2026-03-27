#include "Incompressible/SIMPLE.hpp"
#include "linalg.hpp"
#include "cblas.h"

SIMPLE::SIMPLE(const Domain& domain, const GradientScheme& Ugrads, const GradientScheme& Pgrads, 
    const IterativeSolver& mom_solver, const IterativeSolver& pres_solver, const ConvectionScheme& convScheme, const ViscousModel& viscScheme) 
    : domain(domain), Nc(domain.Nc), Nf(domain.Nf), NBC(domain.NBC), Ugrads(Ugrads), Pgrads(Pgrads), 
    mom_solver(mom_solver), pres_solver(pres_solver), convScheme(convScheme), viscScheme(viscScheme)
{
    iter = 0;

    BCs = new BoundaryCondition*[NBC];

    U[0]       = new double[NDIMS*Nc]();
    G[0]       = new double[NDIMS*NDIMS*Nc]();
    P          = new double[Nc]();
    gradP      = new double[NDIMS*Nc]();

    A_mom      = new double[Nc*Nc]();
    D_mom      = new double[Nc]();
    b_mom[0]   = new double[NDIMS*Nc]();
    rhs_mom[0] = new double[NDIMS*Nc]();
    F          = new double[Nf]();

    A_pres     = new double[Nc*Nc]();
    H_a[0]     = new double[NDIMS*Nc]();
    rhs_pres   = new double[Nf]();
    Pp         = new double[Nc]();

    for (int var = 1; var < NDIMS; var++) {
        U[var]       = U[var-1]       + Nc;
        b_mom[var]   = b_mom[var-1]   + Nc;
        rhs_mom[var] = rhs_mom[var-1] + Nc;
        H_a[var]     = H_a[var-1]     + Nc;

        G[var]       = G[var-1]       + Nc*NDIMS;
    }

    for (int dim = 0; dim < NDIMS; dim++) {
        writer_arrs[dim] = U[dim];
    }
    writer_arrs[NDIMS] = P;
    #if   NDIMS == 2
    std::string varnames[] = {"Velocity_X", "Velocity_Y", "Pressure"};
    #elif NDIMS == 3
    std::string varnames[] = {"Velocity_X", "Velocity_Y", "Velocity_Z", "Pressure"};
    #endif
    writer = new Writer(NDIMS+1, varnames, Nc, domain.cells_permuted, writer_arrs, domain.cperm);
}

SIMPLE::~SIMPLE() {
    for (int bc = 0; bc < NBC; bc++) delete BCs[bc];
    delete[] BCs;

    delete[] U[0];
    delete[] G[0];
    delete[] P;
    delete[] gradP;

    delete[] A_mom;
    delete[] D_mom;
    delete[] b_mom[0];
    delete[] rhs_mom[0];
    delete[] F;

    delete[] A_pres;
    delete[] H_a[0];
    delete[] rhs_pres;
    delete[] Pp;

    delete writer;
}

void SIMPLE::gather_BCs(const int* BCmap, const double UBC[][NDIMS], const double PBC[]) {
    int fcolBC = domain.Ncol-NBC;
    for (int fcol = fcolBC; fcol < domain.Ncol; fcol++) {
        int F1 = domain.mesh.coloff[fcol];
        int F2 = domain.mesh.coloff[fcol+1];
        int bc = fcol-fcolBC;

        int BCtype = BCmap[ domain.mesh.coltag[fcol] ];
        switch (BCtype) {
            case 0: BCs[bc] = new Symmetry(domain); break;
            case 1: BCs[bc] = new Wall(domain)    ; break;
            case 2: BCs[bc] = new Inlet(domain)   ; break;
            case 3: BCs[bc] = new Outlet(domain)  ; break;
        }
        BCs[bc]->init(NDIMS, F1, F2, UBC[bc], PBC[bc], U, P);
    }
}

void SIMPLE::initialise(double P_ref, const double UIC[NDIMS], double PIC) {
    this->P_ref = P_ref;

    for (int dim = 0; dim < NDIMS; dim++) {
        for (int i = 0; i < Nc; i++) {
            U[dim][i] = UIC[dim];
        }
    }

    for (int i = 0; i < Nc; i++) {
        P[i] = PIC;
    }

    update_BC_cells();
    for (int f = 0; f < Nf; f++) {
        int CL = domain.C[f][0];
        int CR = domain.C[f][1];
        
        double ul[NDIMS], ur[NDIMS];

        for (int dim = 0; dim < NDIMS; dim++) {
            ul[dim] = U[dim][CL];
            ur[dim] = U[dim][CR];
        }

        double uf[NDIMS];
        for (int dim = 0; dim < NDIMS; dim++) {
            uf[dim] = ul[dim]*domain.If[f] + ur[dim]*(1.0-domain.If[f]);
        }

        F[f] = dot(domain.Af[f], uf);
    }
}


void SIMPLE::run_solver(int niters) {
    writer->write_soln();
    clear_arr(NDIMS*Nc, gradP);
    Pgrads.compute_gradient(1, P, gradP);
    for (int iter = 0; iter < niters; iter++) {
        solve_mom_eqn();
        solve_pres_eqn();
        writer->write_soln();
        clear_arr(NDIMS*Nc, gradP);
        Pgrads.compute_gradient(1, P, gradP);
        compute_F();
    }
}

void SIMPLE::compute_F() {
    double* If = domain.If;
    for (int f = 0; f < Nf; f++) {
        int CL = domain.C[f][0];
        int CR = domain.C[f][1];

        double gradPf[NDIMS];
        for (int dim = 0; dim < NDIMS; dim++) {
            gradPf[dim] = If[f]*gradP[dim*Nc+CL] + (1.0-If[f])*gradP[dim*Nc+CR]; 
        } 

        double H_af[NDIMS];
        for (int dim = 0; dim < NDIMS; dim++) {
            H_af[dim] = If[f]*H_a[dim][CL] + (1.0-If[f])*H_a[dim][CR]; 
        } 
        
        double iDf = If[f]/D_mom[CL] + (1.0-If[f])/D_mom[CR];

        double Uf[NDIMS];
        for (int dim = 0; dim < NDIMS; dim++) {
            Uf[dim] = H_af[dim] - iDf*gradPf[dim];
        }

        F[f] = dot(domain.Af[f], Uf);
    }
}

void SIMPLE::underrelax() {
    const double alpha_u = 0.8;
    int Nc = domain.Nc;
    int Nf = domain.Nf;
    
    for (int dim = 0; dim < NDIMS; dim++) {
        for (int i = 0; i < Nc; i++) {
            b_mom[dim][i] += A_mom[i+i*Nc]*U[dim][i]*(1.0-alpha_u)/alpha_u;
        }
    }
    for (int i = 0; i < Nc; i++) {
        A_mom[i+i*Nc] *= 1.0/alpha_u;
    }
    
}   

void SIMPLE::solve_mom_eqn() {
    update_BC_cells();

    clear_arr(NDIMS*NDIMS*Nc, G[0]);
    for (int dim = 0; dim < NDIMS; dim++) {
        Ugrads.compute_gradient(NDIMS, U[dim], G[dim]);
    }

    fastclear(Nc, Nf, A_mom, domain.C);
    clear_arr(Nc*NDIMS, b_mom[0]);
    convScheme.discretise_momentum(A_mom, b_mom, F, G, domain);
    viscScheme.discretise_viscosity(A_mom);
    underrelax();
    for (int dim = 0; dim < NDIMS; dim++) {
        for (int i = 0; i < Nc; i++) {
            rhs_mom[dim][i] = b_mom[dim][i] - gradP[dim*Nc+i];
        }
    }
    update_BC_mom_matrix();

    for (int dim = 0; dim < NDIMS; dim++) {
        int niters = mom_solver.solve(A_mom, U[dim], rhs_mom[dim]);
        printf("mom %d iters: %d\n", dim, niters);
    }
}

void SIMPLE::solve_pres_eqn() {
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, Nc, NDIMS, Nc, -1.0, A_mom, Nc, U[0], Nc, 0.0, H_a[0], Nc); // COL MAJOR
    get_diags(Nc, A_mom, D_mom);
    for (int dim = 0; dim < NDIMS; dim++) {
        for (int i = 0; i < Nc; i++) {
            H_a[dim][i] += D_mom[i]*U[dim][i] + b_mom[dim][i];
        }
    }

    for (int i = 0; i < Nc; i++) {
        D_mom[i] = 1.0/D_mom[i];
    }
 
    for (int dim = 0; dim < NDIMS; dim++) {
        for (int i = 0; i < Nc; i++) {
            H_a[dim][i] *= D_mom[i];
        }
    }

    clear_arr(Nc, rhs_pres);
    for (int f = 0; f < Nf; f++) {
        int CL = domain.C[f][0];
        int CR = domain.C[f][1];

        double H_al[NDIMS], H_ar[NDIMS];
        for (int dim = 0; dim < NDIMS; dim++) {
            H_al[dim] = H_a[dim][CL];
            H_ar[dim] = H_a[dim][CR];
        }

        double H_af[NDIMS];
        for (int dim = 0; dim < NDIMS; dim++) {
            H_af[dim] = domain.If[f]*H_al[dim] + (1-domain.If[f])*H_ar[dim];
        }

        double prod = dot(domain.Af[f], H_af);

        rhs_pres[CL] += prod*domain.iv[CL];
        rhs_pres[CR] -= prod*domain.iv[CR];
    }

    fastclear(Nc, Nf, A_pres, domain.C);
    laplacian_nu(Nf, Nc, A_pres, D_mom, domain);
    update_BC_pressure_matrix();
    cblas_dcopy(Nc, P, 1, Pp, 1);
    int niters = pres_solver.solve(A_pres, Pp, rhs_pres);
    printf("pres  iters: %d\n", niters);
    double alpha = 0.2;
    for (int i = 0; i < Nc; i++) {
        P[i] = (1.0-alpha)*P[i] + alpha*Pp[i];
    }
}


void SIMPLE::update_BC_cells() {
    for (int bc = 0; bc < NBC; bc++) {
        BCs[bc]->update_cells();
    }
}

void SIMPLE::update_BC_mom_matrix() {
    for (int bc = 0; bc < NBC; bc++) {
        BCs[bc]->modify_mom_matrix(Nc, A_mom, rhs_mom);
    }
}

void SIMPLE::update_BC_pressure_matrix() {
    for (int bc = 0; bc < NBC; bc++) {
        BCs[bc]->modify_pressure_matrix(Nc, A_pres, rhs_pres);
    }
}



