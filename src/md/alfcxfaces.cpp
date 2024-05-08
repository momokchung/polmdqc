// Author: Moses KJ Chung
// Year:   2024

#include "alfc.h"
#include "alfcxfaces.h"

namespace polmdqc
{
///////////////////////////////////////////////////
//                                               //
//  alfcxfaces  --  generates the list of faces  //
//                                               //
///////////////////////////////////////////////////

// "alfcxfaces" procedure generates the list of boundary
// faces in the alpha complex

void alfcxfaces(std::vector<Tetrahedron>& tetra, std::vector<Face>& faces)
{
    int face_edge[4][3] = {
        { 2, 1, 0},
        { 4, 3, 0},
        { 5, 3, 1},
        { 5, 4, 2},
    };

    // Each triangle is defined implicitly as the interface
    // between two tetrahedra i and j 

    int i, j, k;
    double coef;

    int ntetra = tetra.size();
    faces.clear();
    faces.reserve(ntetra);

    int e_1, e_2, e_3;

    for(int idx = 0; idx < ntetra; idx++) {
        // "dead" tetrahedron are ignored
        if (tetra[idx].info[1]==0) continue;

        for(int itrig = 0; itrig < 4; itrig++) {

            if (tetra[idx].info[2+itrig]==0) continue;

            int jtetra = tetra[idx].neighbors[itrig];

            i = tetra[idx].vertices[other3[itrig][0]];
            j = tetra[idx].vertices[other3[itrig][1]];
            k = tetra[idx].vertices[other3[itrig][2]];

            e_1 = tetra[idx].info_edge[face_edge[itrig][0]];
            e_2 = tetra[idx].info_edge[face_edge[itrig][1]];
            e_3 = tetra[idx].info_edge[face_edge[itrig][2]];

            if (jtetra==-1) {
                coef = 1.0;
                if (tetra[idx].info[6]==1) coef = 0.5;
                Face f = Face(i, j, k, e_1, e_2, e_3, coef);
                faces.push_back(f);
            }
            else if (jtetra > idx) {
                coef = 1.0;
                if (tetra[idx].info[6]==1 && tetra[jtetra].info[6]==1) {
                    coef = 0.0;
                }
                else if (tetra[idx].info[6]==1 || tetra[jtetra].info[6]==1) {
                    coef = 0.5;
                }
                Face f = Face(i, j, k, e_1, e_2, e_3, coef);
                faces.push_back(f);
            }
        }
    }
}
}
