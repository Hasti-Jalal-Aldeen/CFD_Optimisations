#include "CAoS/BoundaryConditions.hpp"
#include <cmath>
using namespace CAoS;

void Symmetry::update_cells(double U[][NVARS]) {
    for (int f = F1; f < F2; f++) {
        int CL = domain.C[0][f];
        int CR = domain.C[1][f];
        
        double n[NDIMS];
        for (int dim = 0; dim < NDIMS; dim++) {
            n[dim] = domain.nf[dim][f];
        }

        double nor = dot(&U[CL][1], n);
        U[CR][0] = U[CL][0];
        for (int dim = 0; dim < NDIMS; dim++) {
            U[CR][dim+1] = U[CL][dim+1] - 2*nor*n[dim];
        }
        U[CR][NDIMS+1] = U[CL][NDIMS+1];
    }
}
