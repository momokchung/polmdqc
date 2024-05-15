// Author: Moses KJ Chung
// Year:   2024

#pragma once
#include "edge.h"
#include "face.h"
#include "mathConst.h"
#include "precision.h"
#include "tetrahedron.h"
#include "vertex.h"
#include <vector>

namespace polmdqc
{
/////////////////////////////////////////
//                                     //
//  alphavol  --  Alpha volumes class  //
//                                     //
/////////////////////////////////////////

class AlphaVol {
public:
    template <bool compder>
    void alphavol(std::vector<Vertex>& vertices, std::vector<Tetrahedron>& tetra,
        std::vector<Edge>& edges, std::vector<Face>& faces,
        real* ballwsurf, real* ballwvol, real* ballwmean, real* ballwgauss,
        real* dsurf_coord, real* dvol_coord, real* dmean_coord, real* dgauss_coord);

private:
    real eps = 1e-14;

    real teteps = 1e-5;

    real twopi = 2 * pi;

    inline real dist2(std::vector<Vertex>& vertices, int n1, int n2);

    inline void twosph(real ra, real ra2, real rb, real rb2,
        real rab, real rab2, real& surfa, real& surfb,
        real& vola, real& volb, real& r, real& phi, real& l);

    template <bool compder>
    inline void twosphder(real ra, real ra2, real rb, real rb2, real rab, real rab2,
        real& surfa, real& surfb, real& vola, real& volb, real& r, real& phi, real& l,
        real& dsurfa, real& dsurfb, real& dvola, real& dvolb, real& dr, real& dphi, real& dl);

    template <bool compder>
    inline void threesphder(real ra, real rb,real rc, real ra2,
        real rb2, real rc2, real rab, real rac, real rbc,
        real rab2, real rac2, real rbc2, real *angle, real deriv[6][3],
        real& surfa, real& surfb, real& surfc, real& vola, real& volb, real& volc,
        real* dsurfa, real* dsurfb, real* dsurfc, real* dvola, real* dvolb, real* dvolc);

    inline real plane_dist(real ra2, real rb2, real rab2);

    template <bool compder>
    inline void tetdihedder(real r12sq, real r13sq, real r14sq,
        real r23sq, real r24sq, real r34sq, real tetvol, real* angle,
        real* cosine, real* sine, real deriv[6][6]);

    template <bool compder>
    inline void tetdihedder3(real r12sq, real r13sq, real r14sq,
        real r23sq, real r24sq, real r34sq, real* angle,
        real* cosine, real* sine, real deriv[6][3]);

    template <bool compder>
    inline void tet3dihedcos(real r12sq, real r13sq, real r14sq,
        real r23sq, real r24sq,real r34sq, real* cosine, real deriv[3][3]);

    template <bool compder>
    inline void tetvorder(real ra2,real rb2,real rc2,real rd2,
        real rab, real rac, real rad, real rbc, real rbd,
        real rcd, real rab2, real rac2, real rad2,real rbc2,
        real rbd2, real rcd2, real tetvol, real* cos_ang, real* sin_ang,
        real deriv[6][6], real& vola, real& volb, real& volc,
        real& vold, real* dvola, real* dvolb, real* dvolc, real* dvold);

    template <bool compder>
    inline real trig_dradius(real a, real b, real c, real *der_r);

    template <bool compder>
    inline real trig_darea(real a, real b, real c, real *der_S);

    inline real sign(real a, real b, real c);

    template <bool compder>
    inline void threesphgss(real ra, real rb, real rc,
        real ra2, real rb2, real rc2, real rab, real rac, real rbc,
        real rab2, real rac2, real rbc2, real& areaA, real& areaB,
        real& areaC, real darea[3][3]);

    inline real tetra_volume(real r12sq, real r13sq, real r14sq, real r23sq, real r24sq, real r34sq);
};
}
