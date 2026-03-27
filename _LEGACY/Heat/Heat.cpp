#include "Heat/Heat.hpp"
#include "linalg.hpp"
#include <string>

BC::BC() {
    return;
}
void BC::init(int F1, int F2, double UBC, int C[][2], double Af[][NDIMS], double* U) {
    this->F1 = F1;
    this->F2 = F2;

    this->C = C;
    this->Af = Af;

    this->U = U;
    this->UBC = UBC;    
}
void BC::update_cells() {
    return;
}
void BC::modify_matrix(int N, double* A, double* rhs) {
    return;
}

void Dirichlet::update_cells() {
    for (int f = F1; f < F2; f++) {
        int CR = C[f][1];
        U[CR] = UBC;
    }
}
void Dirichlet::modify_matrix(int N, double* A, double* rhs) {
    for (int f = F1; f < F2; f++) {
        int CL = C[f][0];
        int CR = C[f][1];

        A[CR+CL*N] = 0.0;
        A[CR+CR*N] = 1.0;
    }

    for (int f = F1; f < F2; f++) {
        int CR = C[f][1];
        rhs[CR] = UBC;
    }
}

void Neumann::update_cells() {
    for (int f = F1; f < F2; f++) {
        int CL = C[f][0];
        int CR = C[f][1];
        U[CR] = U[CL];
    }
}
void Neumann::modify_matrix(int N, double* A, double* rhs) {
    for (int f = F1; f < F2; f++) {
        int CL = C[f][0];
        int CR = C[f][1];

        A[CR+CR*N] =  1.0;
        A[CR+CL*N] = -1.0;

        rhs[CR] = 0.0;
    }
}

Heat::Heat(const Domain& domain, const IterativeSolver& solver) : 
domain(domain), Nc(domain.Nc), Nf(domain.Nf), NBC(domain.NBC), solver(solver)
{
    BCs = new BC*[NBC];

    T = new double[Nc]();
    A = new double[Nc*Nc]();
    rhs = new double[Nc*Nc]();
    std::string varnames[] = {"Temperature"};
    writer = new Writer(1, varnames, Nc, domain.cells_permuted, &T, domain.cperm);
}
Heat::~Heat() {
    for (int bc = 0; bc < NBC; bc++) delete BCs[bc];
    delete[] BCs;
    
    delete[] T;
    delete[] A;
    delete[] rhs;

    delete[] writer;
}
void Heat::initialise(double TIC) {
    for (int i = 0; i < Nc; i++) {
        T[i] = TIC;
    }
}
void Heat::gather_BCs(const int BCmap[], const double TBC[]) {
    int fcolBC = domain.Ncol-NBC;
    for (int fcol = fcolBC; fcol < domain.Ncol; fcol++) {
        int F1 = domain.mesh.coloff[fcol];
        int F2 = domain.mesh.coloff[fcol+1];
        int bc = fcol-fcolBC;

        int BCtype = BCmap[ domain.mesh.coltag[fcol] ];
        switch (BCtype) {
            case 0: BCs[bc] = new Dirichlet(); break;
            case 1: BCs[bc] = new Neumann()  ; break;
        }
        BCs[bc]->init(F1, F2, TBC[bc], domain.C, domain.Af, T);
    }
}
void Heat::run_solver() {
    update_BC_cells();
    writer->write_soln();
    laplacian(Nf, Nc, A, 1.0, domain);
    update_BC_matrix();
    solver.solve(A, T, rhs);
    writer->write_soln();
}
void Heat::update_BC_cells() {
    for (int bc = 0; bc < NBC; bc++) {
        BCs[bc]->update_cells();
    }
}
void Heat::update_BC_matrix() {
    for (int bc = 0; bc < NBC; bc++) {
        BCs[bc]->modify_matrix(Nc, A, rhs);
    }
}