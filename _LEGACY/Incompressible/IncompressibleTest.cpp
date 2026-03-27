#include "Global.hpp"
#include "Mesh.hpp"
#include "Domain.hpp"
#include "Writer.hpp"
#include "linalg.hpp"
#include "GradientSchemes.hpp"
#include "IterativeSolvers.hpp"
#include "Incompressible/ConvectionSchemes.hpp"
#include "Incompressible/Viscosity.hpp"
#include "Incompressible/SIMPLE.hpp"
#include "Heat/Heat.hpp"
#include <stdio.h>

void printmat(int M, int N, double* A) {
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            printf("%f ", A[i*N+j]);
        }
        printf("\n");
    }
    printf("\n");
}

void printmatT(int M, int N, double* A) {
    for (int j = 0; j < N; j++) {
        for (int i = 0; i < M; i++) {
            printf("%f ", A[i*N+j]);
        }
        printf("\n");
    }
    printf("\n");
}


int main(int argc, char** argv) {
    /*
    const int N = 5;
    const double range = 1;
    double M[N*N];
    double A[N*N];
    double x0[N];
    double x1[N];
    double x2[N];
    double b[N];

    // srand(time(0));
    // double mn = 0.7;
    // for (int j = 0; j < N; j++) {
    //     for (int i = 0; i < N; i++) {
    //         double a = range* (double) rand()/RAND_MAX;
    //         A[i+j*N] = (a > 0.5*mn) ? 0.0 : a;
    //     }
    // }

    // for (int i = 0; i < N; i++) {
    //     x0[i] = range* (double) rand()/RAND_MAX;
    //     A[i+i*N] = mn + (range-mn)* (double) rand()/RAND_MAX;
    // }

    // int row = 3;
    // x0[row] = 2.0;
    // for (int j = 0; j < N; j++) {
    //     A[row+j*N] = 0.0;
    // }
    // A[row+row*N] = 1.0;

    srand(time(0));
    double mn = 0.7;
    for (int j = 0; j < N; j++) {
        x0[j] = -1 + 2 * (double) rand()/RAND_MAX;
        for (int i = 0; i < N; i++) {
            M[i+j*N] = -1 + 2 * (double) rand()/RAND_MAX;
        }
    }
    cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, N, N, N, 1.0, M, N, M, N, 0.0, A, N);
    cblas_dgemv(CblasColMajor, CblasNoTrans, N, N, 1.0, A, N, x0, 1, 0.0, b, 1);
    
    int MAXITER = 100;
    GaussSeidel s1 = GaussSeidel(MAXITER, N, 1e-6);
    ConjugateGradient s2 = ConjugateGradient(MAXITER, N, 1e-6);

    int i1 = s1.solve(A, x1, b);
    int i2 = s2.solve(A, x2, b);

    printmat(N, N, A);
    printmat(N, 1, x0);

    printmat(N, 1, x1);
    printmat(N, 1, x2);

    printf("Gauss-seidel: %d/%d, Conj Grad: %d/%d\n", i1, MAXITER, i2, MAXITER);
    */

    const int nvars = NDIMS;
    const int BCmap[4] = {0,1,2,3};
    const double UIC[] = {10,0};
    const double Pref = 100;
    const double PIC  = 0;
    const double nu = 1e-2;
    const double UBC[3][NDIMS] = { {0,0} , {10, 0}, {0,0} };
    const double P_BC[3] = {PIC,PIC,PIC};

    Mesh mesh("../step.hmsh");
    Domain domain(mesh, 0);
    NoGrads Ugrads(domain);
    GaussGreen Pgrads(domain);

    GaussSeidel mom_solver(300, domain.Nc, 1e-6);
    // GaussSeidel pres_solver(300, domain.Nc, 1e-6);
    ConjugateGradient pres_solver(300, domain.Nc, 1e-6);
    
    UD convScheme = UD();
    DNS viscScheme(nu, domain);
    
    SIMPLE system(domain, Ugrads, Pgrads, mom_solver, pres_solver, convScheme, viscScheme);
    system.gather_BCs(BCmap, UBC, P_BC);
    system.initialise(Pref, UIC, PIC);
    system.run_solver(2);
    
    /*
    const int nvars = NDIMS;
    const int BCmap[4] = {0,1,0,0};
    const double TIC = 10;
    const double T_BC[3] = { 10 , 30, 10 };

    Mesh mesh("../step.hmsh");
    Domain domain(mesh, 0);
    ConjugateGradient solver = ConjugateGradient(1000, domain.Nc, 1e-5);

    Heat system = Heat(domain, solver);
    system.initialise(TIC);
    system.gather_BCs(BCmap, T_BC);
    system.run_solver();
    */
    
    return 0;
}