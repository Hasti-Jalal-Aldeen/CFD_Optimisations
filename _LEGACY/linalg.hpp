#pragma once
#include "Global.hpp"
#include "Domain.hpp"

void csrmv(int N, const double* data, const int* cols, const int* off, const double* x, double* y);
void laplacian(int Nf, int Nc, double* A, double nu, const Domain& domain);
void laplacian_nu(int Nf, int Nc, double* A, const double* nu, const Domain& domain);
void clear_arr(int N, double* arr);
void fastclear(int Nc, int Nf, double* A, int C[][2]);
void get_diags(int N, const double* A, double* D);
void permute_csr_rows(int N, const int* perm, const int* con, const int* off, int* con_p, int* off_p);
void permute_csr_cols(int N, const int* iperm, int* con, int* off);

int check_diag_dominance(int N, const double* A);
