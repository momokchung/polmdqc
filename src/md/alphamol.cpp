// Author: Moses KJ Chung and Patrice Koehl
// Year:   2023

#include "alfcx.h"
#include "alfcxedges.h"
#include "alfcxfaces.h"
#include "alfp.h"
#include "alphamol.h"
#include "alphavol.h"
#include "delcx.cpp"
#include "edge.h"
#include "face.h"
#include "inform.h"
#include "tetrahedron.h"
#include "vertex.h"
#include <iostream>
#include <vector>

namespace polmdqc
{
/////////////////////////////////////////////////////
//                                                 //
//  alphamol  --  compute surface area and volume  //
//                                                 //
/////////////////////////////////////////////////////

// "alphamol" computes volume, surface area, mean, and
// gaussian curvature
//
// literature reference:
//
// P. Koehl, A. Akopyan, and H. Edelsbrunner, "Computing the Volume,
// Surface Area, Mean, and Gaussian Curvatures of Molecules and Their
// Derivatives", Journal of Chemical Information and Modeling,
// 63, 973-985, (2023).
//
// github reference:
//
// https://github.com/pkoehl/AlphaMol

void alphamol(int natoms, AlfAtom* alfatoms, real* surf, real* vol, real* mean, real* gauss,
    real* dsurf, real* dvol, real* dmean, real* dgauss, bool deriv)
{
    std::vector<Vertex> vertices;
    std::vector<Tetrahedron> tetra;
    std::vector<Edge> edges;
    std::vector<Face> faces;

    Delcx delcx;

    clock_t start_s,stop_s;
    real total = 0;
    bool alfprint = verbose and alfmeth==AlfMethod::AlphaMol;

    // initialize Delaunay procedure
    if (alfprint) start_s = clock();
    delcx.init(natoms, alfatoms, vertices, tetra);
    if (alfprint) {
        stop_s = clock();
        printf("\n Init compute time         : %10.6f ms\n", (stop_s-start_s)/double(CLOCKS_PER_SEC)*1000);
        total += (stop_s-start_s)/double(CLOCKS_PER_SEC);
    }

    // compute Delaunay triangulation
    if (alfprint) start_s = clock();
    if (alfsos) delcx.regular3D<true>(vertices, tetra);
    else delcx.regular3D<false>(vertices, tetra);
    if (alfprint) {
        stop_s = clock();
        printf("\n Regular3D compute time    : %10.6f ms\n", (stop_s-start_s)/double(CLOCKS_PER_SEC)*1000);
        total += (stop_s-start_s)/double(CLOCKS_PER_SEC);
    }

    // generate alpha complex (with alpha=0.0)
    if (alfprint) start_s = clock();
    real alpha = 0;
    alfcx(vertices, tetra, alpha);
    if (alfprint) {
        stop_s = clock();
        printf("\n AlphaCx compute time      : %10.6f ms\n", (stop_s-start_s)/double(CLOCKS_PER_SEC)*1000);
        total += (stop_s-start_s)/double(CLOCKS_PER_SEC);
    }

    // Compute surface area and, optionally volume of the union of balls.
    // If requested, compute also their derivatives
    if (alfprint) start_s = clock();
    alfcxedges(tetra, edges);
    if (alfprint) {
        stop_s = clock();
        printf("\n AlphaCxEdges compute time : %10.6f ms\n", (stop_s-start_s)/double(CLOCKS_PER_SEC)*1000);
        total += (stop_s-start_s)/double(CLOCKS_PER_SEC);
    }

    if (alfprint) start_s = clock();
    alfcxfaces(tetra, faces);
    if (alfprint) {
        stop_s = clock();
        printf("\n AlphaCxFaces compute time : %10.6f ms\n", (stop_s-start_s)/double(CLOCKS_PER_SEC)*1000);
        total += (stop_s-start_s)/double(CLOCKS_PER_SEC);
    }

    if (alfprint) start_s = clock();
    if (deriv) alphavol<true>(vertices, tetra, edges, faces, surf, vol, mean, gauss, dsurf, dvol, dmean, dgauss);
    else alphavol<false>(vertices, tetra, edges, faces, surf, vol, mean, gauss, dsurf, dvol, dmean, dgauss);
    if (alfprint) {
        stop_s = clock();
        printf("\n Volumes compute time      : %10.6f ms\n", (stop_s-start_s)/double(CLOCKS_PER_SEC)*1000);
        total += (stop_s-start_s)/double(CLOCKS_PER_SEC);
        printf("\n AlphaMol compute time     : %10.6f ms\n", total*1000);
    }
}
}
