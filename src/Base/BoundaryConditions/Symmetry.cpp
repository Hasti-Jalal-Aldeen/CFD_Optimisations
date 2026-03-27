#include "Base/BoundaryConditions.hpp"
#include <cmath>
using namespace Base;

void Symmetry::update_cells(double* U[NVARS]) {
    for (int f = F1; f < F2; f++) {
        int CL = domain.C[f][0];
        int CR = domain.C[f][1];
        
        const double* n = domain.nf[f];

        double ul[NDIMS];
        for (int dim = 0; dim < NDIMS; dim++) {
            ul[dim] = U[dim+1][CL];
        }
        double nor = dot(ul, n);

        U[0][CR] = U[0][CL];
        for (int dim = 0; dim < NDIMS; dim++) {
            U[dim+1][CR] = ul[dim] - 2*nor*n[dim];
        }
        U[NDIMS+1][CR] = U[NDIMS+1][CL];
    }
}
