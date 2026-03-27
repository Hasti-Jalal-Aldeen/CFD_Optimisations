#include "AoSoA/Domain.hpp"
#include <immintrin.h>
using namespace AoSoA;

Domain::Domain(Mesh& mesh) 
: Nf(mesh.Nf), Nc(mesh.Nc), Ncol(mesh.Ncol), NBC(mesh.NBCcol), mesh(mesh)
{
    colour_faces();

    C  = (int    (*)[2]       [VECLEN]) _mm_malloc(sizeof(int   )*Nfb*2*      VECLEN, ALIGN); setarr(0, 2*Nfb*      VECLEN, 0  , C[0][0]    );
    rf = (double (*)[2][NDIMS][VECLEN]) _mm_malloc(sizeof(double)*Nfb*2*NDIMS*VECLEN, ALIGN); setarr(0, 2*Nfb*NDIMS*VECLEN, 0.0, rf[0][0][0]);
    xf = (double (*)   [NDIMS][VECLEN]) _mm_malloc(sizeof(double)*Nfb*  NDIMS*VECLEN, ALIGN); setarr(0,   Nfb*NDIMS*VECLEN, 0.0, xf[0][0]   );
    nf = (double (*)   [NDIMS][VECLEN]) _mm_malloc(sizeof(double)*Nfb*  NDIMS*VECLEN, ALIGN); setarr(0,   Nfb*NDIMS*VECLEN, 0.0, nf[0][0]   );
    Af = (double (*)          [VECLEN]) _mm_malloc(sizeof(double)*Nfb*        VECLEN, ALIGN); setarr(0,   Nfb*VECLEN      , 0.0, Af[0]      );
    If = (double (*)          [VECLEN]) _mm_malloc(sizeof(double)*Nfb*        VECLEN, ALIGN); setarr(0,   Nfb*VECLEN      , 0.0, If[0]      );

    xc = (double (*)   [NDIMS]        ) _mm_malloc(sizeof(double)*(Nc+1)*NDIMS      , ALIGN); setarr(0, (Nc+1)*NDIMS, 0.0, xc[0]);
    iv = (double*                     ) _mm_malloc(sizeof(double)*(Nc+1)            , ALIGN); setarr(0,  Nc+1       , 0.0, iv   );

    get_face_con();
    compute_xf();
    compute_xc();
    compute_rf();
    compute_nf();
    compute_If();
    compute_iv();
}

Domain::~Domain() {
    delete[] coloff;
    _mm_free(C);
    _mm_free(rf);
    _mm_free(xf);
    _mm_free(nf);
    _mm_free(Af);
    _mm_free(If);
    _mm_free(xc);
    _mm_free(iv);
}

void Domain::colour_faces() { // currently only padding for sake of experiment
    coloff = new int[Ncol+1]();

    for (int col = 0; col < Ncol; col++) {
        int nf = mesh.coloff[col+1]-mesh.coloff[col];
        int nb = ceil(nf, VECLEN);
        coloff[col+1] = coloff[col] + nb;
    }
    Nfb = coloff[Ncol];
    int Nfp = Nfb*VECLEN;

    int*  fperm = new int[Nf ]();
    int* ifperm = new int[Nfp];
    setarr(0, Nfb*VECLEN, -1, ifperm);

    for (int col = 0; col < Ncol; col++) {
        int f1  = mesh.coloff[col];
        int f2  = mesh.coloff[col+1];
        int fp1 = coloff[col]*VECLEN;
        for (int f = f1; f < f2; f++) {
            fperm[f] = fp1 + (f-f1);
        }
    }

    for (int f = 0; f < Nf; f++) {
        ifperm[fperm[f]] = f;
    }

    int ncfc = mesh.cf_off[Nc];
    for (int i = 0; i < ncfc; i++) {
        int f = mesh.cf_con[i];
        mesh.cf_con[i] = fperm[f];
    }

    int* fn_buf = new int[2*Nfp]();
    int* fc_buf = new int[2*Nfp];
    setarr(0, 2*Nfp, -1, fc_buf);
    for (int f = 0; f < Nfp; f++) {
        int fp = ifperm[f];

        fn_buf[2*f+0] = (fp > -1) ? mesh.fn_con[2*fp+0] : 0;
        fn_buf[2*f+1] = (fp > -1) ? mesh.fn_con[2*fp+1] : 0;

        fc_buf[2*f+0] = (fp > -1) ? mesh.fc_con[2*fp+0] : Nc;
        fc_buf[2*f+1] = (fp > -1) ? mesh.fc_con[2*fp+1] : Nc;
    }

    delete[] mesh.fn_con;
    delete[] mesh.fc_con;

    mesh.fn_con = fn_buf;
    mesh.fc_con = fc_buf;

    delete[]  fperm;
    delete[] ifperm;
}

void Domain::get_face_con() {
    for (int fb = 0; fb < Nfb; fb++) {
        for (int v = 0; v < VECLEN; v++) {
            int f = fb*VECLEN+v;
            C[fb][0][v] = mesh.fc_con[2*f+0];
            C[fb][1][v] = mesh.fc_con[2*f+1];
        }
    }
}

void Domain::compute_xf() {
    int    n1[VECLEN]         ALIGNED; 
    int    n2[VECLEN]         ALIGNED;
    double nd1[NDIMS][VECLEN] ALIGNED; 
    double nd2[NDIMS][VECLEN] ALIGNED;

    for (int fb = 0; fb < Nfb; fb++) {
        for (int v = 0; v < VECLEN; v++) {
            int f = fb*VECLEN+v;
            n1[v] = mesh.fn_con[2*f+0];
            n2[v] = mesh.fn_con[2*f+1];
        }
        
        for (int dim = 0; dim < NDIMS; dim++) {
            for (int v = 0; v < VECLEN; v++) {
                nd1[dim][v] = mesh.nodes[NDIMS*n1[v]+dim];
                nd2[dim][v] = mesh.nodes[NDIMS*n2[v]+dim];
            }
        }

        for (int dim = 0; dim < NDIMS; dim++) {
            for (int v = 0; v < VECLEN; v++) {
                xf[fb][dim][v] = 0.5*(nd1[dim][v] + nd2[dim][v]);
            }
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
    int FBBC = coloff[Ncol-NBC];
    int CL[VECLEN] ALIGNED;
    int CR[VECLEN] ALIGNED;
    for (int fb = 0; fb < Nfb; fb++) {
        for (int v = 0; v < VECLEN; v++) {
            CL[v] = C[fb][0][v];
            CR[v] = C[fb][1][v];
        }

        for (int dim = 0; dim < NDIMS; dim++) {
            for (int v = 0; v < VECLEN; v++) {
                rf[fb][0][dim][v] = (CL[v] > -1) ? xf[fb][dim][v] - xc[CL[v]][dim] : 0.0;
                rf[fb][1][dim][v] = (CR[v] > -1) ? xf[fb][dim][v] - xc[CR[v]][dim] : 0.0;
            }
        }
    }
    for (int fb = FBBC; fb < Nfb; fb++) {
        for (int dim = 0; dim < NDIMS; dim++) {
            for (int v = 0; v < VECLEN; v++) {
                rf[fb][1][dim][v] = -rf[fb][0][dim][v];
            }
        }
    }
}

void Domain::compute_nf() {
    int    CL[VECLEN]         ALIGNED;
    int    CR[VECLEN]         ALIGNED;
    int    mask[VECLEN]       ALIGNED;

    int    n1[VECLEN]         ALIGNED; 
    int    n2[VECLEN]         ALIGNED;
    double nd1[NDIMS][VECLEN] ALIGNED; 
    double nd2[NDIMS][VECLEN] ALIGNED;

    double n[NDIMS][VECLEN]   ALIGNED;
    double n_dot_rl[VECLEN]   ALIGNED;
    double sgn[VECLEN]        ALIGNED;

    for (int fb = 0; fb < Nfb; fb++) {
        for (int v = 0; v < VECLEN; v++) {
            int f = fb*VECLEN+v;
            n1[v] = mesh.fn_con[2*f+0];
            n2[v] = mesh.fn_con[2*f+1];
        }
        
        for (int dim = 0; dim < NDIMS; dim++) {
            for (int v = 0; v < VECLEN; v++) {
                nd1[dim][v] = mesh.nodes[NDIMS*n1[v]+dim];
                nd2[dim][v] = mesh.nodes[NDIMS*n2[v]+dim];
            }
        }

        for (int v = 0; v < VECLEN; v++) {
            CL[v] = C[fb][0][v];
            CR[v] = C[fb][1][v];
        }
        for (int v = 0; v < VECLEN; v++) {
            mask[v] = (CL[v] > -1) && (CR[v] > -1);
        }
        
        // 2D ONLY
        for (int v = 0; v < VECLEN; v++) {
            n[0][v] = (mask[v]) ?  (nd2[1][v]-nd1[1][v]) : 0.0;
            n[1][v] = (mask[v]) ? -(nd2[0][v]-nd1[0][v]) : 0.0;
        }
        vabsv(n, Af[fb]);

        for (int dim = 0; dim < NDIMS; dim++) {
            for (int v = 0; v < VECLEN; v++) {
                n[dim][v] = (mask[v]) ? n[dim][v]/Af[fb][v] : 0.0;
            }
        }
        dotv(n, rf[fb][0], n_dot_rl);

        for (int v = 0; v < VECLEN; v++) {
            sgn[v] = (n_dot_rl[v] > 0.0) ? 1 : -1;
        }
        for (int v = 0; v < VECLEN; v++) {
            nf[fb][0][v] = sgn[v]*n[0][v];
            nf[fb][1][v] = sgn[v]*n[1][v];
        }
    }
}

void Domain::compute_If() {
    int    CL[VECLEN]         ALIGNED;
    int    CR[VECLEN]         ALIGNED;
    int    mask[VECLEN]       ALIGNED;
    double n_dot_rl[VECLEN]   ALIGNED;
    double n_dot_rr[VECLEN]   ALIGNED;
    for (int fb = 0; fb < Nfb; fb++) {
        dotv(nf[fb], rf[fb][0], n_dot_rl);
        dotv(nf[fb], rf[fb][1], n_dot_rr);

        for (int v = 0; v < VECLEN; v++) {
            CL[v] = C[fb][0][v];
            CR[v] = C[fb][1][v];
        }
        for (int v = 0; v < VECLEN; v++) {
            mask[v] = (CL[v] > -1) && (CR[v] > -1);
        }

        for (int v = 0; v < VECLEN; v++) {
            If[fb][v] = (mask[v]) ? n_dot_rl[v]/(n_dot_rl[v]-n_dot_rr[v]) : 0.5;
        }
    }
}

void Domain::compute_iv() {
    int    CL[VECLEN]         ALIGNED;
    int    CR[VECLEN]         ALIGNED;
    int    mask[VECLEN]       ALIGNED;

    double nf_dot_rl[VECLEN]  ALIGNED;
    double nf_dot_rr[VECLEN]  ALIGNED;

    for (int fb = 0; fb < Nfb; fb++) {
        dotv(nf[fb], rf[fb][0], nf_dot_rl);
        dotv(nf[fb], rf[fb][1], nf_dot_rr);

        for (int v = 0; v < VECLEN; v++) {
            CL[v] = C[fb][0][v];
            CR[v] = C[fb][1][v];
        }

        for (int v = 0; v < VECLEN; v++) {
            iv[CL[v]] += Af[fb][v]*nf_dot_rl[v]/NDIMS;
            iv[CR[v]] -= Af[fb][v]*nf_dot_rr[v]/NDIMS;
        }
    }

    for (int i = 0; i < Nc; i++) iv[i] = 1.0/iv[i];

    for (int fb = 0; fb < Nfb; fb++) {
        for (int v = 0; v < VECLEN; v++) {
            CL[v] = C[fb][0][v];
            CR[v] = C[fb][1][v];
        }
        for (int v = 0; v < VECLEN; v++) {
            iv[CL[v]] = (CL[v] > -1) ? iv[CL[v]] : 0.0;
            iv[CR[v]] = (CR[v] > -1) ? iv[CR[v]] : 0.0;
        }
    }

    int FBBC = coloff[Ncol-NBC];
    for (int fb = FBBC; fb < Nfb; fb++) {
        for (int v = 0; v < VECLEN; v++) {
            CL[v] = C[fb][0][v];
            CR[v] = C[fb][1][v];
        }
        for (int v = 0; v < VECLEN; v++) {
            iv[CR[v]] = iv[CL[v]];
        }
    }
}

void Domain::zero_bc_iv() {
    int FBC = mesh.coloff[Ncol-NBC];
    for (int fb = FBC; fb < Nfb; fb++) {
        for (int v = 0; v < VECLEN; v++) {
            int CR = C[fb][1][v];
            iv[CR] = 0.0;
        }
    }
}

