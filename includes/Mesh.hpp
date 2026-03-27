#pragma once
#include <string>

class Mesh {
    public:
    int dim, Nn, Nc, Nf, Ncol, NBCcol, N_pnames;
    std::string* pnames;
    int *coltag, *coloff;
    double *nodes;
    int *fn_con, *fc_con;
    int *cf_con, *cf_off;
    int *cn_off, *cn_con;
    int *cc_off, *cc_con;

    Mesh(const char* fname);
    ~Mesh();

    void parse_pnames(int Nchars, const char* pnames_raw);
    void colour_faces();

    private:
    void compute_cf_con();
};