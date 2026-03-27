#include "FAoS_CAoS/Domain.hpp"
using namespace FAoS_CAoS;

Domain::Domain(Mesh& mesh) 
: Nf(mesh.Nf), Nc(mesh.Nc), Ncol(mesh.Ncol), NBC(mesh.NBCcol), mesh(mesh)
{
    C  = new int[Nf][2]();
    rf = new double[Nf][2][NDIMS]();
    xf = new double[Nf][NDIMS]();
    nf = new double[Nf][NDIMS]();
    Af = new double[Nf]();
    If = new double[Nf]();

    xc = new double[Nc][NDIMS]();
    iv = new double[Nc]();

    get_face_con();
    compute_xf();
    compute_xc();
    compute_rf();
    compute_nf();
    compute_If();
    compute_iv();
}

Domain::~Domain() {
    delete[] C;
    delete[] rf;
    delete[] xf;
    delete[] nf;
    delete[] Af;
    delete[] If;
    delete[] xc;
    delete[] iv;
}

void Domain::get_face_con() {
    for (int f = 0; f < Nf; f++) {
        C[f][0] = mesh.fc_con[2*f+0];
        C[f][1] = mesh.fc_con[2*f+1];
    }
}

void Domain::compute_xf() {
    for (int f = 0; f < Nf; f++) {
        int n1 = mesh.fn_con[2*f+0];
        int n2 = mesh.fn_con[2*f+1];

        double nd1[NDIMS], nd2[NDIMS];
        for (int dim = 0; dim < NDIMS; dim++) {
            nd1[dim] = mesh.nodes[NDIMS*n1+dim];
            nd2[dim] = mesh.nodes[NDIMS*n2+dim];
        }

        for (int dim = 0; dim < NDIMS; dim++) {
            xf[f][dim] = 0.5*(nd1[dim] + nd2[dim]);
        }
    }
}

void Domain::compute_xc() {
    for (int i = 0; i < Nc; i++) {
        int nnds = mesh.cn_off[i+1]-mesh.cn_off[i];
        for (int n = mesh.cn_off[i]; n < mesh.cn_off[i+1]; n++) {
            int nd = mesh.cn_con[n];
            for (int dim = 0; dim < NDIMS; dim++) {
                xc[i][dim] += mesh.nodes[nd*NDIMS+dim]/nnds;
            }
        }
    }
}

void Domain::compute_rf() {
    int F_BC = mesh.coloff[Ncol-NBC];
    for (int f = 0; f < Nf; f++) {
        int CL = C[f][0];
        int CR = C[f][1];

        for (int dim = 0; dim < NDIMS; dim++) {
            rf[f][0][dim] = xf[f][dim] - xc[CL][dim];
            rf[f][1][dim] = xf[f][dim] - xc[CR][dim];
        }
    }
    for (int f = F_BC; f < Nf; f++) {
        for (int dim = 0; dim < NDIMS; dim++) {
            rf[f][1][dim] = -rf[f][0][dim];
        }
    }
}

void Domain::compute_nf() {
    for (int f = 0; f < Nf; f++) {
        int n1 = mesh.fn_con[2*f+0];
        int n2 = mesh.fn_con[2*f+1];

        double nd1[NDIMS], nd2[NDIMS];
        for (int dim = 0; dim < NDIMS; dim++) {
            nd1[dim] = mesh.nodes[NDIMS*n1+dim];
            nd2[dim] = mesh.nodes[NDIMS*n2+dim];
        }
        
        // 2D ONLY
        nf[f][0] =  (nd2[1]-nd1[1]);
        nf[f][1] = -(nd2[0]-nd1[0]);
        Af[f]    = vabs(nf[f]);
        for (int dim = 0; dim < NDIMS; dim++) {
            nf[f][dim] /= Af[f];
        }

        double nf_dot_rl = dot(nf[f], rf[f][0]);

        int sgn  = (nf_dot_rl > 0) ? 1 : -1;
        nf[f][0] = sgn*nf[f][0];
        nf[f][1] = sgn*nf[f][1];
    }
}

void Domain::compute_If() {
    for (int f = 0; f < Nf; f++) {
        double nf_dot_rL =  dot(nf[f], rf[f][0]);
        double nf_dot_rR = -dot(nf[f], rf[f][1]);

        If[f] = nf_dot_rL/(nf_dot_rL+nf_dot_rR);
    }
}

void Domain::compute_iv() {
    for (int f = 0; f < Nf; f++) {
        double Af_dot_rL = Af[f]*dot(nf[f], rf[f][0]);
        double Af_dot_rR = Af[f]*dot(nf[f], rf[f][1]);

        int CL = C[f][0];
        int CR = C[f][1];
        iv[CL] += Af_dot_rL/NDIMS;
        iv[CR] -= Af_dot_rR/NDIMS;
    }

    for (int i = 0; i < Nc; i++) iv[i] = 1.0/iv[i];

    int F_BC = mesh.coloff[Ncol-NBC];
    for (int f = F_BC; f < Nf; f++) {
        int CL = C[f][0];
        int CR = C[f][1];
        iv[CR] = iv[CL];
    }
}

void Domain::zero_bc_iv() {
    int F_BC = mesh.coloff[Ncol-NBC];
    for (int f = F_BC; f < Nf; f++) {
        int CR = C[f][1];
        iv[CR] = 0.0;
    }
}

