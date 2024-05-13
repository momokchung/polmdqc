// Author: Moses KJ Chung
// Year:   2024

#pragma once
#include "edge.h"
#include "face.h"
#include "tetrahedron.h"
#include "vertex.h"
#include <vector>

namespace polmdqc
{
//////////////////////////////////////
//                                  //
//  alfcx  --  Alpha complex class  //
//                                  //
//////////////////////////////////////

class Alfcx {
public:
    void alfcx(std::vector<Vertex>& vertices, std::vector<Tetrahedron>& tetra,
        std::vector<Edge>& edges, std::vector<Face>& faces, real alpha);

private:
    void alfcxedges(std::vector<Tetrahedron>& tetra, std::vector<Edge>& edges);

    void alfcxfaces(std::vector<Tetrahedron>& tetra, std::vector<Face>& faces);

    inline int findedge(Tetrahedron t, int i1, int j1);

    void alf_tetra(real* a, real* b, real* c, real* d,
        real ra, real rb, real rc, real rd, int& iflag, real alpha);

    void alf_trig(real* a, real* b, real* c, real* d, real* e,
        real ra,real rb, real rc, real rd, real re, int ie,
        int& irad,int& iattach, real alpha);

    void alf_edge(std::vector<Vertex>& vertices, real* a, real* b, real ra, real rb,
        real* cg, std::vector<int>& listcheck, int& irad, int& iattach, real alpha);

    inline void edge_radius(real* a, real* b, real ra, real rb,
        real* Dab, real* Sab, real* Tab, int& testr, real alpha);

    inline void edge_attach(real* a, real* b, real* c, real ra, real rb,
        real rc, real* Dab, real* Sab, real* Tab, int& testa);

    inline void triangle_attach(real* a, real* b, real* c, real* d,
        real ra, real rb, real rc, real rd, real S[3][4], real T[2][3],
        real Dabc, int& testa, int& memory);

    inline void triangle_radius(real* a, real* b, real* c, real ra, real rb,
        real rc, real S[3][4], real T[2][3], real Dabc, int& testr, real alpha, int& memory);

    inline void vertex_attach(real* a, real* b, real ra, real rb, int& testa, int& testb);

    inline void get_coord2(std::vector<Vertex>& vertices, int ia, int ja, real* a, real* b, real* cg, real& ra, real& rb);

    inline void get_coord4(std::vector<Vertex>& vertices, int ia, int ja, int ka, int la, real* a, real* b,
        real* c, real* d, real& ra, real& rb, real& rc, real& rd);

    inline void get_coord5(std::vector<Vertex>& vertices, int ia, int ja, int ka, int la,
        int ma, real* a, real* b, real* c, real* d, real* e,
        real& ra, real& rb, real& rc, real& rd, real& re);

    // ALFCX_GMP alf_gmp;

protected:
    int other3[4][3] = {
        {1, 2, 3},
        {0, 2, 3},
        {0, 1, 3},
        {0, 1, 2} };

    int face_edge[4][3] = {
        { 2, 1, 0},
        { 4, 3, 0},
        { 5, 3, 1},
        { 5, 4, 2},
    };

    int face_info[6][2] = {
        {0, 1},
        {0, 2},
        {0, 3},
        {1, 2},
        {1, 3},
        {2, 3}};

    int face_pos[6][2] = {
        {1, 0},
        {2, 0},
        {3, 0},
        {2, 1},
        {3, 1},
        {3, 2}};

    int pair[6][2] = {
        {2, 3},
        {1, 3},
        {1, 2},
        {0, 3},
        {0, 2},
        {0, 1}};

    real eps;
};
}
