// Author: Moses KJ Chung and Patrice Koehl
// Year:   2023

#include "alfcx.h"
#include "alfcxedges.h"
#include "alfcxfaces.h"
#include "alfp.h"
#include "alphamol.h"
#include "alphavol.h"
#include "delaunay.h"
#include "edge.h"
#include "face.h"
#include "inform.h"
#include "initdelcx.h"
#include "tetrahedron.h"
#include "vertex.h"
#include <iostream>
#include <queue>
#include <stack>
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

void alphamol(int natoms, AlfAtom* alfatoms, real& wsurf, real& wvol, real& wmean, real& wgauss,
    real* surf, real* vol, real* mean, real* gauss,
    real* dsurf, real* dvol, real* dmean, real* dgauss, bool deriv)
{
    std::vector<Vertex> vertices;
    std::vector<Tetrahedron> tetra;
    std::vector<Edge> edges;
    std::vector<Face> faces;
    std::queue<std::pair<int,int>> link_facet;
    std::queue<std::pair<int,int>> link_index;
    std::stack<int> free;
    std::vector<int> kill;

    clock_t start_s,stop_s;
    real total = 0;

    // initialize Delaunay procedure
    start_s = clock();
    initdelcx(natoms, alfatoms, vertices, tetra, link_facet, link_index, free, kill);
    stop_s = clock();

    if (verbose and alfmeth==AlfMethod::AlphaMol) {
        printf("\n Initdelcx compute time    : %10.6f ms\n", (stop_s-start_s)/double(CLOCKS_PER_SEC)*1000);
        total += (stop_s-start_s)/double(CLOCKS_PER_SEC);
    }

    // compute Delaunay triangulation
    start_s = clock();
    delaunay(vertices, tetra, link_facet, link_index, free, kill);
    stop_s = clock();

    if (verbose and alfmeth==AlfMethod::AlphaMol) {
        printf("\n Delaunay compute time     : %10.6f ms\n", (stop_s-start_s)/double(CLOCKS_PER_SEC)*1000);
        total += (stop_s-start_s)/double(CLOCKS_PER_SEC);
    }

    // generate alpha complex (with alpha=0.0)
    start_s = clock();
    real alpha = 0;
    alfcx(vertices, tetra, alpha);
    stop_s = clock();

    if (verbose and alfmeth==AlfMethod::AlphaMol) {
        printf("\n AlphaCx compute time      : %10.6f ms\n", (stop_s-start_s)/double(CLOCKS_PER_SEC)*1000);
        total += (stop_s-start_s)/double(CLOCKS_PER_SEC);
    }

    // Compute surface area and, optionally volume of the union of balls.
    // If requested, compute also their derivatives
    start_s = clock();
    alfcxedges(tetra, edges);
    stop_s = clock();

    if (verbose and alfmeth==AlfMethod::AlphaMol) {
        printf("\n AlphaCxEdges compute time : %10.6f ms\n", (stop_s-start_s)/double(CLOCKS_PER_SEC)*1000);
        total += (stop_s-start_s)/double(CLOCKS_PER_SEC);
    }

    start_s = clock();
    alfcxfaces(tetra, faces);
    stop_s = clock();

    if (verbose and alfmeth==AlfMethod::AlphaMol) {
        printf("\n AlphaCxFaces compute time : %10.6f ms\n", (stop_s-start_s)/double(CLOCKS_PER_SEC)*1000);
        total += (stop_s-start_s)/double(CLOCKS_PER_SEC);
    }

    start_s = clock();
    alphavol(vertices, tetra, edges, faces, wsurf, wvol, wmean, wgauss, surf, vol, mean, gauss,
        dsurf, dvol, dmean, dgauss, deriv);
    stop_s = clock();

    if (verbose and alfmeth==AlfMethod::AlphaMol) {
        printf("\n Volumes compute time      : %10.6f ms\n", (stop_s-start_s)/double(CLOCKS_PER_SEC)*1000);
        total += (stop_s-start_s)/double(CLOCKS_PER_SEC);
    }

    if (verbose and alfmeth==AlfMethod::AlphaMol) {
        printf("\n AlphaMol compute time     : %10.6f ms\n", total*1000);
    }
}
}
