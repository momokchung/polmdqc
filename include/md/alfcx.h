/* ====================================================================
   alfcx
 
	builds the alpha complex based on the weighted Delaunay triangulation

 ==================================================================== */

#pragma once

#include "alfcxgmp.h"
#include "edge.h"
#include "face.h"
#include "tetrahedron.h"
#include <iostream>
#include <vector>

/* ====================================================================
  class
 ==================================================================== */

class ALFCX {
public:

    void alfcx(double alpha, std::vector<Vertex>& vertices, std::vector<Tetrahedron>& tetra);

    void alphacxEdges(std::vector<Tetrahedron>& tetra, std::vector<Edge>& edges);

    void alphacxFaces(std::vector<Tetrahedron>& tetra, std::vector<Face>& faces);


private:

    int findEdge(Tetrahedron t, int i1, int j1);

    void alf_tetra(double *a, double *b, double *c, double *d,
    double ra, double rb, double rc, double rd, int *iflag, double alpha);

    void alf_trig(double *a, double *b, double *c, double *d, double *e,
    double ra,double rb, double rc, double rd, double re, int ie, 
    int *irad,int *iattach, double alpha);

    void alf_edge(std::vector<Vertex>& vertices, double *a, double *b, double ra, 
    double rb, double *cg, std::vector<int>& listcheck, int *irad, int *iattach, 
    double alpha);

    void edge_radius(double *a, double *b, double ra, double rb,
    double *Dab, double *Sab, double *Tab, int *testr, double alpha,
    int *memory);

    void edge_attach(double *a, double *b, double *c, double ra, double rb,
    double rc, double *Dab, double *Sab, double *Tab, int *testa, int *memory);

    void triangle_attach(double *a, double *b, double *c, double *d,
    double ra, double rb, double rc, double rd, double S[3][4], double T[2][3],
    double Dabc, int *testa, int *memory);

    void triangle_radius(double *a, double *b, double *c, double ra, double rb,
    double rc, double S[3][4], double T[2][3], double Dabc, int *testr, double alpha,
    int *memory);

    void vertex_attach(double *a, double *b, double ra, double rb, int *testa,
    int *testb);

    void get_coord2(std::vector<Vertex>& vertices, int ia, int ja,
        double *a, double *b, double *cg, double *ra, double *rb);

    void get_coord4(std::vector<Vertex>& vertices, int ia, int ja, int ka, int la,
        double *a, double *b, double *c, double *d, double *ra, 
        double *rb, double *rc, double *rd);

    void get_coord5(std::vector<Vertex>& vertices, int ia, int ja, int ka, int la,
        int ma, double *a, double *b, double *c, double *d, double *e,
        double *ra, double *rb, double *rc, double *rd, double *re);

protected:

    int other3[4][3] = {
        {1, 2, 3},
        {0, 2, 3},
        {0, 1, 3},
        {0, 1, 2} };

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

    double eps = 1.e-5;
};
