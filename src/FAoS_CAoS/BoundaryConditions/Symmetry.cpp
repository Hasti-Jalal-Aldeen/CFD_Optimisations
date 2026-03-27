#include "FAoS_CAoS/BoundaryConditions.hpp"
#include <cmath>
using namespace FAoS_CAoS;

void Symmetry::update_cells(double U[][NVARS]) {
    for (int f = F1; f < F2; f++) {
        int CL = domain.C[f][0];
        int CR = domain.C[f][1];

        double nor = dot(&U[CL][1], domain.nf[f]);
        U[CR][0] = U[CL][0];
        for (int dim = 0; dim < NDIMS; dim++) {
            U[CR][dim+1] = U[CL][dim+1] - 2*nor*domain.nf[f][dim];
        }
        U[CR][NDIMS+1] = U[CL][NDIMS+1];
    }
}
