// Author: Moses KJ Chung
// Year:   2024

#pragma once
#include "libfunc.h"
#include "mathConst.h"
#include "precision.h"
#include "vertex.h"
#include "tetra.h"
#include <vector>

namespace polmdqc
{
////////////////////////////////////////////////////////////////
//                                                            //
//  alfsphere  --  compute surf & vol of overlapping spheres  //
//                                                            //
////////////////////////////////////////////////////////////////

constexpr real aseps = 1e-14;

inline real dist2(std::vector<Vertex>& vertices, int n1, int n2)
{
    real x;
    real dist = 0;
    for(int i = 0; i < 3; i++) {
        x = vertices[n1].coord[i] - vertices[n2].coord[i];
        dist += x*x;
    }

    return dist;
}

/////////////////////////////////////////////////////////
//                                                     //
//  twosph  --  surf & vol of two overlapping spheres  //
//                                                     //
/////////////////////////////////////////////////////////

// "twosph" calculates the surface area and volume of the
// intersection of two spheres

inline void twosph(real ra, real ra2, real rb, real rb2,
    real rab, real rab2, real& surfa, real& surfb,
    real& vola, real& volb, real& r, real& phi, real& l)
{
    real cosine,vala,valb,lambda,ha,hb;
    real Aab,sa,ca,sb,cb;
    constexpr real twopi = 2 * pi;

    // Get distance between center of sphere A and Voronoi plane
    // between A and B
    lambda = plane_dist(ra2, rb2, rab2);
    valb = lambda*rab;
    vala = rab-valb;

    // Get height of the cap of sphere A occluded by sphere B
    ha = ra - vala;

    // same for sphere B ...
    hb = rb - valb;

    // get surfaces of intersection
    surfa = twopi*ra*ha;
    surfb = twopi*rb*hb;

    // now get volume
    Aab = pi*(ra2-vala*vala);

    sa = ra*(surfa);
    ca = vala*Aab;

    vola = (sa-ca)/3;

    sb = rb*(surfb);
    cb = valb*Aab;

    volb = (sb-cb)/3;

    // get radius of the circle of intersection between the two spheres
    r = REAL_SQRT(ra2 - vala*vala);

    // get angle between normals of the sphere at a point on this circle
    cosine = (ra2+rb2-rab2)/(2.0*ra*rb);
    if (REAL_ABS(cosine - 1) < aseps) cosine = 1;
    else if (REAL_ABS(cosine + 1) < aseps) cosine = -1;
    phi = REAL_ACOS(cosine);

    l = vala/ra + valb/rb;
}

///////////////////////////////////////////////////////////////////
//                                                               //
//  twosphder  --  surf & vol derivs of two overlapping spheres  //
//                                                               //
///////////////////////////////////////////////////////////////////

// "twosphder" calculates the surface area and volume derivatives
// of the intersection of two spheres

inline void twosphder(real ra, real ra2, real rb, real rb2, real rab, real rab2,
    real& surfa, real& surfb, real& vola, real& volb, real& r, real& phi, real& l,
    real& dsurfa, real& dsurfb, real& dvola, real& dvolb, real& dr, real& dphi, real& dl, bool compder)
{
    real cosine,vala,valb,lambda,ha,hb;
    real Aab,sa,ca,sb,cb;
    real dera,derb;
    constexpr real twopi = 2 * pi;

    // Get distance between center of sphere A and Voronoi plane
    // between A and B
    lambda = plane_dist(ra2, rb2, rab2);
    valb = lambda*rab;
    vala = rab-valb;

    // get height of the cap of sphere A occluded by sphere B
    ha = ra - vala;

    // same for sphere B ...
    hb = rb - valb;

    // get surfaces of intersection
    surfa = twopi*ra*ha;
    surfb = twopi*rb*hb;

    // now get volume
    Aab = pi*(ra2-vala*vala);

    sa = ra*(surfa);
    ca = vala*Aab;

    vola = (sa-ca)/3;

    sb = rb*(surfb);
    cb = valb*Aab;

    volb = (sb-cb)/3;

    // get radius of the circle of intersection between the two spheres
    r = REAL_SQRT(ra2 - vala*vala);

    // get angle between normals of the sphere at a point on this circle
    cosine = (ra2+rb2-rab2)/(2.0*ra*rb);
    if (REAL_ABS(cosine - 1) < aseps) cosine = 1;
    else if (REAL_ABS(cosine + 1) < aseps) cosine = -1;
    phi = REAL_ACOS(cosine);
    l = vala/ra + valb/rb;

    if (!compder) return;

    dera = - lambda;
    derb = lambda - 1;

    dsurfa = twopi*ra*dera;
    dsurfb = twopi*rb*derb;

    dvola = -Aab*lambda;
    dvolb = -(dvola) - Aab;

    dr   = -vala*lambda/(r);
    dphi = rab/(ra*rb*REAL_SQRT(1-cosine*cosine));
    dl   = lambda/ra + (1.0-lambda)/rb;
}

///////////////////////////////////////////////////////////////////////
//                                                                   //
//  threesphder  --  surf & vol derivs of three overlapping spheres  //
//                                                                   //
///////////////////////////////////////////////////////////////////////

// "threesphder" calculates the surface area and volume derivatives
// of the intersection of three spheres

inline void threesphder(real ra, real rb,real rc, real ra2,
    real rb2, real rc2, real rab, real rac, real rbc,
    real rab2, real rac2, real rbc2, real *angle, real deriv[6][3],
    real& surfa, real& surfb, real& surfc, real& vola, real& volb, real& volc,
    real* dsurfa, real* dsurfb, real* dsurfc, real* dvola, real* dvolb, real* dvolc, bool compder)
{
    real a1,a2,a3,s2,c1,c2;
    real seg_ang_ab,seg_ang_ac,seg_ang_bc;
    real ang_dih_ap,ang_dih_bp,ang_dih_cp;
    real val1,val2,val3,l1,l2,l3;
    real val1b,val2b,val3b;
    real ang_abc,ang_acb,ang_bca;
    real cos_abc,cos_acb,cos_bca;
    real sin_abc,sin_acb,sin_bca;
    real s_abc,s_acb,s_bca;
    real rho_ab2,rho_ac2,rho_bc2;
    real drho_ab2,drho_ac2,drho_bc2;
    real val_abc,val_acb,val_bca;
    real val2_abc,val2_acb,val2_bca;
    real der_val1b,der_val1,der_val2b,der_val2,der_val3b,der_val3;
    real cosine[6],sine[6];
    constexpr real twopi = 2 * pi;

    l1 = plane_dist(ra2, rb2, rab2);
    l2 = plane_dist(ra2, rc2, rac2);
    l3 = plane_dist(rb2, rc2, rbc2);

    val1 = l1*rab; val2 = l2*rac; val3 = l3*rbc;
    val1b = rab - val1; val2b = rac - val2; val3b = rbc - val3;

    // We consider the tetrahedron (A,B,C,P) where P is the
    // point of intersection of the three spheres such that (A,B,C,P) is ccw.
    // The edge lengths in this tetrahedron are: rab, rac, rAP=ra, rbc, rBP=rb, rCP=rc

    tetdihedder3(rab2, rac2, ra2, rbc2, rb2, rc2, angle, cosine, sine, deriv, compder);

    // the seg_ang_ are the dihedral angles around the three edges AB, AC and BC

    seg_ang_ab = angle[0];
    seg_ang_ac = angle[1];
    seg_ang_bc = angle[3];

    // the ang_dih_ are the dihedral angles around the three edges AP, BP and CP
    ang_dih_ap = angle[2];
    ang_dih_bp = angle[4];
    ang_dih_cp = angle[5];

    a1 = ra*(1-2*ang_dih_ap);
    a2 = 2*seg_ang_ab*val1b;
    a3 = 2*seg_ang_ac*val2b;

    surfa = twopi*ra*(a1 - a2 - a3);

    a1 = rb*(1-2*ang_dih_bp);
    a2 = 2*seg_ang_ab*val1;
    a3 = 2*seg_ang_bc*val3b;

    surfb = twopi*rb*(a1 - a2 - a3);

    a1 = rc*(1-2*ang_dih_cp);
    a2 = 2*seg_ang_ac*val2;
    a3 = 2*seg_ang_bc*val3;

    surfc = twopi*rc*(a1 - a2 - a3);

    // compute volumes of the three caps
    ang_abc = twopi*seg_ang_ab;
    ang_acb = twopi*seg_ang_ac;
    ang_bca = twopi*seg_ang_bc;

    cos_abc = cosine[0];
    sin_abc = sine[0];
    cos_acb = cosine[1];
    sin_acb = sine[1];
    cos_bca = cosine[3];
    sin_bca = sine[3];

    rho_ab2 = ra2 - val1b*val1b;
    rho_ac2 = ra2 - val2b*val2b;
    rho_bc2 = rb2 - val3b*val3b;

    val_abc = ang_abc - sin_abc*cos_abc; s_abc = rho_ab2*val_abc;
    val_acb = ang_acb - sin_acb*cos_acb; s_acb = rho_ac2*val_acb;
    val_bca = ang_bca - sin_bca*cos_bca; s_bca = rho_bc2*val_bca;

    s2 = ra*(surfa);
    c1 = val1b*s_abc;
    c2 = val2b*s_acb;

    vola = (s2 - c1 - c2)/3;

    s2 = rb*(surfb);
    c1 = val1*s_abc;
    c2 = val3b*s_bca;

    volb = (s2 - c1 - c2)/3;

    s2 = rc*(surfc);
    c1 = val2*s_acb;
    c2 = val3*s_bca;

    volc = (s2 - c1 - c2)/3;

    if (!compder) return;

    der_val1b = l1; der_val1  = 1-l1;
    der_val2b = l2; der_val2  = 1-l2;
    der_val3b = l3; der_val3  = 1-l3;

    dsurfa[0] = -2*ra*(
        twopi*seg_ang_ab*der_val1b +
        (ra*deriv[2][0] +
        val1b*deriv[0][0] +val2b*deriv[1][0]));
    dsurfa[1] = -2*ra*(
        twopi*seg_ang_ac*der_val2b +
        (ra*deriv[2][1] +
        val1b*deriv[0][1] +val2b*deriv[1][1]));
    dsurfa[2] = -2*ra*( ra*deriv[2][2] +
        val1b*deriv[0][2]+val2b*deriv[1][2]);

    dsurfb[0] = -2*rb*(
        twopi*seg_ang_ab*der_val1
        +(rb*deriv[4][0]+
        val1*deriv[0][0]+val3b*deriv[3][0]));
    dsurfb[1] = -2*rb*(rb*deriv[4][1]+
        val1*deriv[0][1]+val3b*deriv[3][1]);
    dsurfb[2] = -2*rb*(
        twopi*seg_ang_bc*der_val3b
        +(rb*deriv[4][2]+
        val1*deriv[0][2]+val3b*deriv[3][2]));

    dsurfc[0] = -2*rc*(rc*deriv[5][0]+
            val2*deriv[1][0]+val3*deriv[3][0]);
    dsurfc[1] = -2*rc*(
        twopi*seg_ang_ac*der_val2
        +(rc*deriv[5][1]+
        val2*deriv[1][1]+val3*deriv[3][1]));
    dsurfc[2] = -2*rc*(
        twopi*seg_ang_bc*der_val3
        +(rc*deriv[5][2]+
        val2*deriv[1][2]+val3*deriv[3][2]));

    drho_ab2 = -2*der_val1b*val1b;
    drho_ac2 = -2*der_val2b*val2b;
    drho_bc2 = -2*der_val3b*val3b;

    val2_abc = rho_ab2*(1 - cos_abc*cos_abc + sin_abc*sin_abc);
    val2_acb = rho_ac2*(1 - cos_acb*cos_acb + sin_acb*sin_acb);
    val2_bca = rho_bc2*(1 - cos_bca*cos_bca + sin_bca*sin_bca);

    dvola[0] = ra*dsurfa[0] - der_val1b*s_abc - 
        (val1b*deriv[0][0]*val2_abc + val2b*deriv[1][0]*val2_acb)
        - val1b*drho_ab2*val_abc;
    dvola[0] = dvola[0]/3;
    dvola[1] = ra*dsurfa[1] - der_val2b*s_acb - 
        (val1b*deriv[0][1]*val2_abc + val2b*deriv[1][1]*val2_acb)
        - val2b*drho_ac2*val_acb;
    dvola[1] = dvola[1]/3;
    dvola[2] = ra*dsurfa[2] - 
        (val1b*deriv[0][2]*val2_abc + val2b*deriv[1][2]*val2_acb);
    dvola[2] = dvola[2]/3;

    dvolb[0] = rb*dsurfb[0] - der_val1*s_abc - 
        (val1*deriv[0][0]*val2_abc + val3b*deriv[3][0]*val2_bca)
        - val1*drho_ab2*val_abc;
    dvolb[0] = dvolb[0]/3;
    dvolb[1] = rb*dsurfb[1] - 
        (val1*deriv[0][1]*val2_abc + val3b*deriv[3][1]*val2_bca);
    dvolb[1] = dvolb[1]/3;
    dvolb[2] = rb*dsurfb[2] - der_val3b*s_bca - 
        (val1*deriv[0][2]*val2_abc + val3b*deriv[3][2]*val2_bca)
        - val3b*drho_bc2*val_bca;
    dvolb[2] = dvolb[2]/3;

    dvolc[0] = rc*dsurfc[0] - 
        (val2*deriv[1][0]*val2_acb + val3*deriv[3][0]*val2_bca);
    dvolc[0] = dvolc[0]/3;
    dvolc[1] = rc*dsurfc[1] - der_val2*s_acb - 
        (val2*deriv[1][1]*val2_acb + val3*deriv[3][1]*val2_bca)
        - val2*drho_ac2*val_acb;
    dvolc[1] = dvolc[1]/3;
    dvolc[2] = rc*dsurfc[2] - der_val3*s_bca - 
        (val2*deriv[1][2]*val2_acb + val3*deriv[3][2]*val2_bca)
        - val3*drho_bc2*val_bca;
    dvolc[2] = dvolc[2]/3;
}
}
