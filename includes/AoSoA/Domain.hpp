#pragma once
#include "Mesh.hpp"
#include "Global.hpp"
namespace AoSoA {

class Domain {
    public:
    int Nfb, Ncol;
    const int Nf, Nc, NBC;

    // Mesh
    Mesh& mesh;
    int* coloff;

    // Faces
    int    (*C )[2][VECLEN];
    double (*rf)[2][NDIMS][VECLEN]; // r[face][left/right][dim]
    double (*nf)[NDIMS][VECLEN];
    double (*xf)[NDIMS][VECLEN];
    double (*Af)[VECLEN];
    double (*If)[VECLEN];

    // Cells
    double (*xc)[NDIMS];
    double  *iv;

    Domain(Mesh& mesh);
    ~Domain();
    void zero_bc_iv();

    private:
    void colour_faces();
    void compute_xc();
    void get_face_con();
    void compute_xf();
    void compute_rf();
    void compute_nf();
    void compute_If();
    void compute_iv();
};

}