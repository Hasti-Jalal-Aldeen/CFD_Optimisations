#include "Incompressible/Viscosity.hpp"
#include "linalg.hpp"

ViscousModel::ViscousModel(double nu, const Domain& domain) : nu(nu), domain(domain) {
    return;
}

void ViscousModel::discretise_viscosity(double* A) const {
    return;
}

void DNS::discretise_viscosity(double* A) const {
    laplacian(domain.Nf, domain.Nc, A, -nu, domain);
}