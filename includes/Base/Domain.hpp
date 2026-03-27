#pragma once
#include "Mesh.hpp"
#include "Global.hpp"
namespace Base {

class Domain {
    public:
    const int Nf,Nc;
    const int Ncol, NBC;

    // Mesh
    Mesh& mesh;

    // Faces
    int    (*C )[2];
    double (*rf)[2][NDIMS]; // r[face][left/right][dim]
    double (*nf)[NDIMS];
    double (*xf)[NDIMS];
    double  *Af;
    double  *If;

    // Cells
    double (*xc)[NDIMS];
    double  *iv;

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