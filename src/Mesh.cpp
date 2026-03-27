#include "Mesh.hpp"
#include "Global.hpp"
#include <cstdio>
#include <iostream>
#include <string>

Mesh::Mesh(const char* fname) {
    FILE* file;
    file = fopen(fname, "rb");
    
    int Ndat[8], Nchars;
    fread(Ndat, sizeof(int), 8, file);
    dim      = Ndat[0];
    Nn       = Ndat[1];
    Ncol     = Ndat[2];
    NBCcol   = Ndat[3];
    Nc       = Ndat[4];
    Nf       = Ndat[5];
    N_pnames = max(1,Ndat[6]);
    Nchars   = Ndat[7];
    
    char* pnames_raw = new char[Nchars];
    pnames  = new std::string[N_pnames];
    coltag  = new int[Ncol];
    coloff  = new int[Ncol+1];
    nodes   = new double[Nn*dim];
    fn_con  = new int[Nf*2];
    fc_con  = new int[Nf*2];
    cn_off  = new int[Nc+1];
    cc_off  = new int[Nc+1];
    
    fread(pnames_raw, sizeof(char)  , Nchars  , file);
    fread(coltag    , sizeof(int)   , Ncol    , file);
    fread(coloff    , sizeof(int)   , Ncol+1  , file);
    fread(nodes     , sizeof(double), Nn*dim  , file);
    fread(fn_con    , sizeof(int)   , 2*Nf    , file);
    fread(fc_con    , sizeof(int)   , 2*Nf    , file);
    fread(cn_off    , sizeof(int)   , Nc+1    , file);
    int N_cn_con = cn_off[Nc];
    cn_con = new int[N_cn_con];
    fread(cn_con    , sizeof(int)   , N_cn_con, file);
    fread(cc_off    , sizeof(int)   , Nc+1    , file);
    int N_cc_con = cc_off[Nc];
    cc_con = new int[N_cc_con];
    fread(cc_con    , sizeof(int)   , N_cc_con, file);
    fclose(file);

    parse_pnames(Nchars, pnames_raw);

    delete[] pnames_raw;
    compute_cf_con();
}

Mesh::~Mesh() {
    delete[] pnames;
    delete[] coltag;
    delete[] coloff;
    delete[] nodes;
    delete[] cn_con;
    delete[] cn_off;
    delete[] fn_con;
    delete[] fc_con;
    delete[] cc_off;
    delete[] cc_con;

    delete[] cf_off;
    delete[] cf_con;
}

void Mesh::parse_pnames(int Nchars, const char* pnames_raw) {
    int off[N_pnames];
    off[0] = 0;

    for (int str = 0; str < N_pnames-1; str++) {
        for (int chr = off[str]; chr < Nchars; chr++) {
            if (pnames_raw[chr] == '\0') {
                off[str+1] = chr+1;
                break;
            }
        }
    }

    for (int str = 0; str < N_pnames; str++) {
        pnames[str] = &pnames_raw[off[str]];
    }
}

void Mesh::compute_cf_con() {
    cf_off = new int[Nc+1]();
    for (int f = 0; f < Nf; f++) {
        int CL = fc_con[2*f+0];
        int CR = fc_con[2*f+1];

        cf_off[CL+1]++;
        cf_off[CR+1]++;
    }
    for (int i = 0; i < Nc; i++) {
        cf_off[i+1] += cf_off[i];
    }

    cf_con = new int[cf_off[Nc]]();
    int* counter = new int[Nc]();

    for (int f = 0; f < Nf; f++) {
        int CL = fc_con[2*f+0];
        int CR = fc_con[2*f+1];

        int posL = cf_off[CL] + counter[CL];
        counter[CL]++;
        int posR = cf_off[CR] + counter[CR];
        counter[CR]++;
        
        cf_con[posL] = f;
        cf_con[posR] = f;
    }
    delete[] counter;
}

void Mesh::colour_faces() {
    
}
