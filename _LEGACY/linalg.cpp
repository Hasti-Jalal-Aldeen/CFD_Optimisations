#include "linalg.hpp"
#include "Global.hpp"
#include <math.h>
#include <cblas.h>

void csrmv(int N, const double* data, const int* cols, const int* off, const double* x, double* y) {
    for (int i = 0; i < N; i++) {
        y[i] = 0.0;

        int colptr1 = off[i];
        int colptr2 = off[i+1];

        for (int colptr = colptr1; colptr < colptr2; colptr++) {
            int col = cols[colptr];
            y[i] += data[colptr]*x[col];
        }
    }
}

void laplacian(int Nf, int Nc, double* A, double nu, const Domain& domain) {
    for (int f = 0; f < Nf; f++) {
        int CL = domain.C[f][0];
        int CR = domain.C[f][1];

        double d[NDIMS];
        for (int dim = 0; dim < NDIMS; dim++) {
            d[dim] = domain.rf[f][0][dim] - domain.rf[f][1][dim];
        }

        double A_d  = domain.Af[f] /sqrt( dot(d, d) );

        A[CL+CL*Nc] -= nu*A_d*domain.iv[CL];
        A[CL+CR*Nc] += nu*A_d*domain.iv[CL];
        A[CR+CR*Nc] -= nu*A_d*domain.iv[CR];
        A[CR+CL*Nc] += nu*A_d*domain.iv[CR];
    }
}
// DIVIDE BY VOLUME
void laplacian_nu(int Nf, int Nc, double* A, const double* nu, const Domain& domain) {
    for (int f = 0; f < Nf; f++) {
        int CL = domain.C[f][0];
        int CR = domain.C[f][1];

        double d[NDIMS];
        for (int dim = 0; dim < NDIMS; dim++) {
            d[dim] = domain.rf[f][0][dim] - domain.rf[f][1][dim];
        }

        double A_d  = domain.Af[f] /sqrt( dot(d, d) );
        double nu_f = domain.If[f]*nu[CL] + (1-domain.If[f])*nu[CR];

        A[CL+CL*Nc] -= nu_f*A_d*domain.iv[CL];
        A[CL+CR*Nc] += nu_f*A_d*domain.iv[CL];
        A[CR+CR*Nc] -= nu_f*A_d*domain.iv[CR];
        A[CR+CL*Nc] += nu_f*A_d*domain.iv[CR];
    }
}

void clear_arr(int N, double* arr) {
    for (int i = 0; i < N; i++) {
        arr[i] = 0.0;
    }
}

void fastclear(int Nc, int Nf, double* A, int C[][2]) {
    // clear diagonals
    for (int i = 0; i < Nc; i++) {
        A[i+i*Nc] = 0.0;
    }

    // clear off-diagonals
    for (int f = 0; f < Nf; f++) {
        int CL = C[f][0];
        int CR = C[f][1];

        A[CL+CR*Nc] = 0.0;
        A[CR+CL*Nc] = 0.0;
    }
}

void get_diags(int N, const double* A, double* D) {
    for (int i = 0; i < N; i++) {
        D[i] = A[i+i*N];
    }
}

void permute_csr_rows(int N, const int* perm, const int* con, const int* off, int* con_p, int* off_p) {
    off_p[0] = 0;
    for (int row_n = 0; row_n < N; row_n++) {
        int row_o = perm[row_n];
        int ncols = off[row_o+1]-off[row_o];
        off_p[row_n+1] = off_p[row_n]+ncols;

        int col_n = off_p[row_n];
        for (int col_o = off[row_o]; col_o < off[row_o+1]; col_o++) {
            con_p[col_n] = con[col_o];
            col_n++;
        }
    }
}

void permute_csr_cols(int N, const int* iperm, int* con, int* off) {
    for (int i = 0; i < N; i++) {
        for (int j = off[i]; j < off[i+1]; j++) {
            con[j] = iperm[con[j]];
        }

        for (int j = off[i]+1; j < off[i+1]; j++) {
            int key = con[j];
            int k = j-1;

            while ( (k >= off[i]) && (con[k] > key) ) {
                con[k+1] = con[k];
                k--;
            }
            con[k+1] = key;
        }
    }
}

int check_diag_dominance(int N, const double* A) {
    for (int i = 0; i < N; i++) {
        double off_diags = 0.0;
        for (int j = 0; j < N; j++) {
            off_diags += fabs(A[i+j*N]);
        }
        double abs_aii = fabs(A[i+i*N]);
        off_diags -= abs_aii;
        if (abs_aii < off_diags) return i;
    }
    return -1;
} 


