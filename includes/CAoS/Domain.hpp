#pragma once
#include "Mesh.hpp"
#include "Global.hpp"
namespace CAoS {

class Domain {
    public:
    const int Nf, Nfp, Nc;
    const int Ncol, NBC;

    // Mesh
    Mesh& mesh;

    // Faces
    int*    C[2];
    double* rf[NDIMS][2]; // r[dim][left/right][face]
    double* nf[NDIMS];
    double* xf[NDIMS];
    double* Af;
    double* If;

    // Cells
    double* xc[NDIMS];
    double* iv;

    Domain(Mesh& mesh);
    ~Domain();
    void zero_bc_iv();

    private:
    void compute_xc();
    void get_face_con();
    void compute_xf();
    void compute_rf();
    void compute_nf();
    void compute_If();
    void compute_iv();
};

}