#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "Global.hpp"
#include "Mesh.hpp"

#include "Base/Domain.hpp"
#include "Base/GradientSchemes.hpp"
#include "Base/Limiters.hpp"
#include "Base/Flux.hpp"
#include "Base/BoundaryConditions.hpp"
#include "Base/TimeIntegrators.hpp"

#include "SoA/Domain.hpp"
#include "SoA/GradientSchemes.hpp"
#include "SoA/Limiters.hpp"
#include "SoA/Flux.hpp"
#include "SoA/BoundaryConditions.hpp"
#include "SoA/TimeIntegrators.hpp"

#include "VEC/Domain.hpp"
#include "VEC/GradientSchemes.hpp"
#include "VEC/Limiters.hpp"
#include "VEC/Flux.hpp"
#include "VEC/BoundaryConditions.hpp"
#include "VEC/TimeIntegrators.hpp"

#include "CAoS/Domain.hpp"
#include "CAoS/GradientSchemes.hpp"
#include "CAoS/Limiters.hpp"
#include "CAoS/Flux.hpp"
#include "CAoS/BoundaryConditions.hpp"
#include "CAoS/TimeIntegrators.hpp"

#include "FAoS_CAoS/Domain.hpp"
#include "FAoS_CAoS/GradientSchemes.hpp"
#include "FAoS_CAoS/Limiters.hpp"
#include "FAoS_CAoS/Flux.hpp"
#include "FAoS_CAoS/BoundaryConditions.hpp"
#include "FAoS_CAoS/TimeIntegrators.hpp"

#include "AoSoA/Domain.hpp"
#include "AoSoA/GradientSchemes.hpp"
#include "AoSoA/Limiters.hpp"
#include "AoSoA/Flux.hpp"
#include "AoSoA/BoundaryConditions.hpp"
#include "AoSoA/TimeIntegrators.hpp"

int main(int argc, char** argv) {
    const int nvars = 4;
    const double gma = 1.4;

    double M = 2.0;
    double p = 1.0/(gma*M*M);
    const double QIC[] = {1,1,0,p};
    const double QBC[][NVARS] = { {1,1,0,p} , {1,1,0,p}, {1,1,0,p}, {1,1,0,p}};
    double dt = 1e-3;

    Mesh mesh("../wedge.hmsh");

    int niter = 1000;

    const Base::BCtype BCmap[] = {Base::SYMMETRY, Base::SYMMETRY, Base::SYMMETRY, Base::SUPINLET, Base::SUPOUTLET};
    Base::Domain domain(mesh);
    Base::GaussGreen grads(domain);
    Base::Limiter limiter(domain);
    Base::Flux flux(nvars, gma, domain);
    Base::Euler intg(nvars, domain, flux, grads, limiter);
    intg.gather_BCs(BCmap, QBC);
    intg.initialise(QIC);
    clock_t t1 = clock();
    intg.run(niter*dt, dt, niter*dt);
    clock_t t2 = clock();
    
    const SoA::BCtype BCmapSoA[] = {SoA::SYMMETRY, SoA::SYMMETRY, SoA::SYMMETRY, SoA::SUPINLET, SoA::SUPOUTLET};
    SoA::Domain domainSoA(mesh);
    SoA::GaussGreen gradsSoA(domainSoA);
    SoA::Limiter limiterSoA(domainSoA);
    SoA::Flux fluxSoA(nvars, gma, domainSoA);
    SoA::Euler intgSoA(nvars, domainSoA, fluxSoA, gradsSoA, limiterSoA);
    intgSoA.gather_BCs(BCmapSoA, QBC);
    intgSoA.initialise(QIC);
    clock_t t3 = clock();
    intgSoA.run(niter*dt, dt, niter*dt);
    clock_t t4 = clock();
    
    const VEC::BCtype BCmapVEC[] = {VEC::SYMMETRY, VEC::SYMMETRY, VEC::SYMMETRY, VEC::SUPINLET, VEC::SUPOUTLET};
    VEC::Domain domainVEC(mesh);
    VEC::GaussGreen gradsVEC(domainVEC);
    VEC::Limiter limiterVEC(domainVEC);
    VEC::Flux fluxVEC(nvars, gma, domainVEC);
    VEC::Euler intgVEC(nvars, domainVEC, fluxVEC, gradsVEC, limiterVEC);
    intgVEC.gather_BCs(BCmapVEC, QBC);
    intgVEC.initialise(QIC);
    clock_t t5 = clock();
    intgVEC.run(niter*dt, dt, niter*dt);
    clock_t t6 = clock();
    
    const CAoS::BCtype BCmapCAoS[] = {CAoS::SYMMETRY, CAoS::SYMMETRY, CAoS::SYMMETRY, CAoS::SUPINLET, CAoS::SUPOUTLET};
    CAoS::Domain domainCAoS(mesh);
    CAoS::GaussGreen gradsCAoS(domainCAoS);
    CAoS::Limiter limiterCAoS(domainCAoS);
    CAoS::Flux fluxCAoS(nvars, gma, domainCAoS);
    CAoS::Euler intgCAoS(nvars, domainCAoS, fluxCAoS, gradsCAoS, limiterCAoS);
    intgCAoS.gather_BCs(BCmapCAoS, QBC);
    intgCAoS.initialise(QIC);
    clock_t t7 = clock();
    intgCAoS.run(niter*dt, dt, niter*dt);
    clock_t t8 = clock();

    const FAoS_CAoS::BCtype BCmapFAoSCAoS[] = {FAoS_CAoS::SYMMETRY, FAoS_CAoS::SYMMETRY, FAoS_CAoS::SYMMETRY, FAoS_CAoS::SUPINLET, FAoS_CAoS::SUPOUTLET};
    FAoS_CAoS::Domain domainFAoSCAoS(mesh);
    FAoS_CAoS::GaussGreen gradsFAoSCAoS(domainFAoSCAoS);
    FAoS_CAoS::Limiter limiterFAoSCAoS(domainFAoSCAoS);
    FAoS_CAoS::Flux fluxFAoSCAoS(nvars, gma, domainFAoSCAoS);
    FAoS_CAoS::Euler intgFAoSCAoS(nvars, domainFAoSCAoS, fluxFAoSCAoS, gradsFAoSCAoS, limiterFAoSCAoS);
    intgFAoSCAoS.gather_BCs(BCmapFAoSCAoS, QBC);
    intgFAoSCAoS.initialise(QIC);
    clock_t t9 = clock();
    intgFAoSCAoS.run(niter*dt, dt, niter*dt);
    clock_t t10 = clock();

    const AoSoA::BCtype BCmapAoSoA[] = {AoSoA::SYMMETRY, AoSoA::SYMMETRY, AoSoA::SYMMETRY, AoSoA::SUPINLET, AoSoA::SUPOUTLET};
    AoSoA::Domain domainAoSoA(mesh);
    AoSoA::GaussGreen gradsAoSoA(domainAoSoA);
    AoSoA::Limiter limiterAoSoA(domainAoSoA);
    AoSoA::Flux fluxAoSoA(nvars, gma, domainAoSoA);
    AoSoA::Euler intgAoSoA(nvars, domainAoSoA, fluxAoSoA, gradsAoSoA, limiterAoSoA);
    intgAoSoA.gather_BCs(BCmapAoSoA, QBC);
    intgAoSoA.initialise(QIC);
    clock_t t11 = clock();
    intgAoSoA.run(niter*dt, dt, niter*dt);
    clock_t t12 = clock();

    double f = 2000;
    double st1 = f * (double) (t2-t1)/CLOCKS_PER_SEC;
    double st2 = f * (double) (t4-t3)/CLOCKS_PER_SEC;
    double st3 = f * (double) (t6-t5)/CLOCKS_PER_SEC;
    double st4 = f * (double) (t8-t7)/CLOCKS_PER_SEC;
    double st5 = f * (double) (t10-t9)/CLOCKS_PER_SEC;
    double st6 = f * (double) (t12-t11)/CLOCKS_PER_SEC;

    printf("f AoS: %f f SoA: %f f VEC: %f c AoS: %f f AoS c AoS %f AoSoA %f\n", st1, st2, st3, st4, st5, st6);

    return 0;
}
