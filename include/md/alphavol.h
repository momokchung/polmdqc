// Author: Moses KJ Chung
// Year:   2024

#pragma once
#include "edge.h"
#include "face.h"
#include "precision.h"
#include "tetrahedron.h"
#include "vertex.h"
#include <vector>

namespace polmdqc
{
///////////////////////////////////////////////////////////////
//                                                           //
//  alphavol  --  compute surf, vol, mean and gaussian curv  //
//                                                           //
///////////////////////////////////////////////////////////////

template <bool compder>
void alphavol(std::vector<Vertex>& vertices, std::vector<Tetrahedron>& tetra,
    std::vector<Edge>& edges, std::vector<Face>& faces,
    real* ballwsurf, real* ballwvol, real* ballwmean, real* ballwgauss,
    real* dsurf_coord, real* dvol_coord, real* dmean_coord, real* dgauss_coord);
}
