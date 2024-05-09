// Author: Moses KJ Chung
// Year:   2024

#pragma once
#include "alfsphere.h"
#include "libfunc.h"
#include "mathConst.h"

namespace polmdqc
{
//////////////////////////////////////////////////////
//                                                  //
//  tetra  --  geometric properties of tetrahedron  //
//                                                  //
//////////////////////////////////////////////////////

constexpr real teteps = 1e-14;

inline real plane_dist(real ra2, real rb2, real rab2)
{
    real lambda = 0.50 - (ra2-rb2)/(2*rab2);

    return lambda;
}

////////////////////////////////////////////////////////////////
//                                                            //
//  tetdihedder  --  tetrahedron dihedral angles derivatives  //
//                                                            //
////////////////////////////////////////////////////////////////

// "tetdihedder" computes the derivative of the six
// dihedral angles of atetrahedronfrom its edge lengths

template <bool compder>
inline void tetdihedder(real r12sq, real r13sq, real r14sq,
    real r23sq, real r24sq, real r34sq, real* angle,
    real* cosine, real* sine, real deriv[6][6])
{
    real val1,val2,val3,val4,vala;
    real val123,val124,val134,val234;
    real val213,val214,val314,val324,val312;
    real det12,det13,det14,det23,det24,det34;
    real minori[4];
    real dminori[4][6] = {0};
    real det[6],dnum[6][6],val[4];
    real dist[6];
    constexpr real twopi = 2 * pi;

    // Define the Cayley Menger matrix:
    // M = ( 0      r12^2  r13^2  r14^2  1)
    //     ( r12^2  0      r23^2  r24^2  1)
    //     ( r13^2  r23^2  0      r34^2  1)
    //     ( r14^2  r24^2  r34^2  0      1)
    //     ( 1      1      1      1      0)
    // Compute all minors M(i,i): determinant of the Cayley-Menger matrix with row i
    // and column j removed
    // These determinants are of the form:
    // det = | 0 a b 1 |
    //       | a 0 c 1 |
    //       | b c 0 1 |
    //       | 1 1 1 0 |
    // then:
    // det = (c - a - b )^2 - 4ab

    val234 = (r34sq - r23sq - r24sq);
    val134 = (r34sq - r14sq - r13sq);
    val124 = (r24sq - r12sq - r14sq);
    val123 = (r23sq - r12sq - r13sq);

    minori[0] = val234*val234 - 4*r23sq*r24sq;
    minori[1] = val134*val134 - 4*r13sq*r14sq;
    minori[2] = val124*val124 - 4*r12sq*r14sq;
    minori[3] = val123*val123 - 4*r12sq*r13sq;

    val4 = 1.0/REAL_SQRT(-minori[0]);
    val3 = 1.0/REAL_SQRT(-minori[1]);
    val2 = 1.0/REAL_SQRT(-minori[2]);
    val1 = 1.0/REAL_SQRT(-minori[3]);

    if constexpr (compder) val[0] = val4; val[1] = val3; val[2] = val2; val[3] = val1;

    // Now compute all angles (in fact, cosine of the angle):
    //           (-1)^(i+j) * det(Mij) 
    // cos(i,j)= ---------------------
    //            sqrt(M(i,i)*M(j,j))
    // where det(Mij) = M(i,j) is the determinant of the Cayley-Menger matrix with row i
    // and column j removed

    det12 = -2*r12sq*val134 - val123*val124;
    det13 = -2*r13sq*val124 - val123*val134;
    det14 = -2*r14sq*val123 - val124*val134;

    val213 = r13sq -r12sq -r23sq;
    val214 = r14sq -r12sq -r24sq;
    val312 = r12sq -r13sq -r23sq;
    val314 = r14sq -r13sq -r34sq;
    val324 = r24sq -r23sq -r34sq;

    det23 = -2*r23sq*val214 - val213*val234;
    det24 = -2*r24sq*val213 - val214*val234;
    det34 = -2*r34sq*val312 - val314*val324;

    cosine[0] = det12*val1*val2;
    cosine[1] = det13*val1*val3;
    cosine[2] = det14*val2*val3;
    cosine[3] = det23*val1*val4;
    cosine[4] = det24*val2*val4;
    cosine[5] = det34*val3*val4;

    for (int i = 0; i < 6; i++) {
        if (REAL_ABS(cosine[i] - 1) < teteps) cosine[i] = 1;
        else if (REAL_ABS(cosine[i] + 1) < teteps) cosine[i] = -1;
    }

    for (int i = 0; i < 6; i++) {
        angle[i] = REAL_ACOS(cosine[i]);
        sine[i]  = REAL_SIN(angle[i]);
        angle[i] /= twopi;
    }

    if constexpr (!compder) return;

    det[5] = det12; det[4] = det13; det[3] = det14;
    det[2] = det23; det[1] = det24; det[0] = det34;

    dist[0] = REAL_SQRT(r12sq); dist[1] = REAL_SQRT(r13sq); dist[2] = REAL_SQRT(r14sq);
    dist[3] = REAL_SQRT(r23sq); dist[4] = REAL_SQRT(r24sq); dist[5] = REAL_SQRT(r34sq);

    // Now compute derivatives of the angles with respect to the edge lengths
    // Since (see above):
    //                       num(i,j)
    // cos(ang(i,j)) = --------------------
    //                  sqrt(M(i,i)*M(j,j))

    //   d(ang(i,j))                       dnum(i,j)                             (M(i,i)dM(j,j) +M(j,j)*dM(i,i))
    // ------------sin(ang(i,j)) =  -----------------------   -0.5*num(i,j) -----------------------------------------
    //    dr(a,b)                   sqrt(M(i,i)M(j,j))dr(a,b)                    M(i,i)M(j,j) sqrt(M(i,i)M(j,j))

    // which we can rewrite as:

    //   d(ang(i,j))                cosine(i,j)   dnum(i,j)                    dM(j,j) +  dM(i,i))
    // ------------sin(ang(i,j)) = -----------  -----------  -0.5*cosine(i,j)( -------- + ---------)
    //    dr(a,b)                    num(i,j)     dr(a,b)                       M(j,j)      M(i,i)

    dminori[0][3] = -(val234 + 2*r24sq); dminori[0][4] = -(val234 + 2*r23sq); dminori[0][5] = val234;
    dminori[1][1] = -(val134 + 2*r14sq); dminori[1][2] = -(val134 + 2*r13sq); dminori[1][5] = val134;
    dminori[2][0] = -(val124 + 2*r14sq); dminori[2][2] = -(val124 + 2*r12sq); dminori[2][4] = val124;
    dminori[3][0] = -(val123 + 2*r13sq); dminori[3][1] = -(val123 + 2*r12sq); dminori[3][3] = val123;

    dnum[5][0] = -2*val134+val123+val124; dnum[5][1] = 2*r12sq + val124; dnum[5][2] = 2*r12sq + val123;
    dnum[5][3] = -val124; dnum[5][4] = -val123; dnum[5][5] = -2*r12sq;

    dnum[4][0] = 2*r13sq+val134; dnum[4][1] = -2*val124 + val123 + val134; dnum[4][2] = 2*r13sq + val123;
    dnum[4][3] = -val134; dnum[4][4] = -2*r13sq; dnum[4][5] = -val123;

    dnum[3][0] = 2*r14sq+val134; dnum[3][1] = 2*r14sq + val124; dnum[3][2] = -2*val123 + val124 + val134;
    dnum[3][3] = -2*r14sq; dnum[3][4] = -val134; dnum[3][5] = -val124;

    dnum[2][0] = 2*r23sq+val234; dnum[2][1] = -val234; dnum[2][2] = -2*r23sq;
    dnum[2][3] = -2*val214+val213+val234; dnum[2][4] = 2*r23sq+val213; dnum[2][5] = -val213;

    dnum[1][0] = 2*r24sq+val234; dnum[1][1] = -2*r24sq; dnum[1][2] = -val234;
    dnum[1][3] = 2*r24sq+val214; dnum[1][4] = -2*val213 + val214 + val234; dnum[1][5] = -val214;

    dnum[0][0] = -2*r34sq; dnum[0][1] = 2*r34sq+val324; dnum[0][2] = -val324;
    dnum[0][3] = 2*r34sq+val314; dnum[0][4] = -val314; dnum[0][5] = -2*val312 + val314 + val324;

    int k = 0;
    int jj;
    for (int i = 0; i < 3; i++) {
        for (int j = i+1; j < 4; j++) {
            jj = 5-k;
            if (det[k] != 0) {
                vala = cosine[jj]/sine[jj];
                val1 = -vala/det[k];
                val2 = vala/minori[j];
                val3 = vala/minori[i];
                for (int l = 0; l < 6; l++) {
                    deriv[jj][l] = val1*dnum[k][l]+val2*dminori[j][l]+val3*dminori[i][l];
                    deriv[jj][l] *= 2*dist[l];
                }
            }
            else {
                vala = -val[i]*val[j]/sine[jj];
                for (int l = 0; l < 6; l++) {
                    deriv[jj][l] = vala*dnum[k][l];
                    deriv[jj][l] *= 2*dist[l];
                }
            }
            k++;
        }
    }
}

/////////////////////////////////////////////////////////////////
//                                                             //
//  tetdihedder3  --  compute dihedral angles and derivatives  //
//                                                             //
/////////////////////////////////////////////////////////////////

// "tetdihedder3" computes the six dihedral angles of a tetrahedron
// ABCD from its edge lengths as well as their derivatives with
// respect to the 3 edge lengths AB, AC and BC

template <bool compder>
inline void tetdihedder3(real r12sq, real r13sq, real r14sq,
    real r23sq, real r24sq, real r34sq, real* angle,
    real* cosine, real* sine, real deriv[6][3])
{
    real val1,val2,val3,val4,vala;
    real val123,val124,val134,val234;
    real val213,val214,val314,val324,val312;
    real det12,det13,det14,det23,det24,det34;
    real minori[4];
    real dminori[4][3] = {0};
    real dist[3],det[6],dnum[6][3],val[4];
    constexpr real twopi = 2 * pi;

    // Define the Cayley Menger matrix:
    // M = ( 0      r12^2  r13^2  r14^2  1)
    //     ( r12^2  0      r23^2  r24^2  1)
    //     ( r13^2  r23^2  0      r34^2  1)
    //     ( r14^2  r24^2  r34^2  0      1)
    //     ( 1      1      1      1      0)
    // Compute all minors M(i,i): determinant of the Cayley-Menger matrix with row i
    // and column j removed
    // These determinants are of the form:
    // det = | 0 a b 1 |
    //       | a 0 c 1 |
    //       | b c 0 1 |
    //       | 1 1 1 0 |
    // then:
    // det = (c - a - b )^2 - 4ab

    val234 = (r34sq - r23sq - r24sq);
    val134 = (r34sq - r14sq - r13sq);
    val124 = (r24sq - r12sq - r14sq);
    val123 = (r23sq - r12sq - r13sq);

    minori[0] = val234*val234 - 4*r23sq*r24sq;
    minori[1] = val134*val134 - 4*r13sq*r14sq;
    minori[2] = val124*val124 - 4*r12sq*r14sq;
    minori[3] = val123*val123 - 4*r12sq*r13sq;

    val4 = 1.0/REAL_SQRT(-minori[0]);
    val3 = 1.0/REAL_SQRT(-minori[1]);
    val2 = 1.0/REAL_SQRT(-minori[2]);
    val1 = 1.0/REAL_SQRT(-minori[3]);

    val[0] = val4; val[1] = val3; val[2] = val2; val[3] = val1;

    // Now compute all angles (in fact, cosine of the angle):
    //           (-1)^(i+j) * det(Mij) 
    // cos(i,j)= ---------------------
    //            sqrt(M(i,i)*M(j,j))
    // where det(Mij) = M(i,j) is the determinant of the Cayley-Menger matrix with row i
    // and column j removed

    det12 = -2*r12sq*val134 - val123*val124;
    det13 = -2*r13sq*val124 - val123*val134;
    det14 = -2*r14sq*val123 - val124*val134;

    val213 = r13sq -r12sq -r23sq;
    val214 = r14sq -r12sq -r24sq;
    val312 = r12sq -r13sq -r23sq;
    val314 = r14sq -r13sq -r34sq;
    val324 = r24sq -r23sq -r34sq;

    det23 = -2*r23sq*val214 - val213*val234;
    det24 = -2*r24sq*val213 - val214*val234;
    det34 = -2*r34sq*val312 - val314*val324;

    cosine[0] = det12*val1*val2;
    cosine[1] = det13*val1*val3;
    cosine[2] = det14*val2*val3;
    cosine[3] = det23*val1*val4;
    cosine[4] = det24*val2*val4;
    cosine[5] = det34*val3*val4;

    for (int i = 0; i < 6; i++) {
        if (REAL_ABS(cosine[i] - 1) < teteps) cosine[i] = 1;
        else if (REAL_ABS(cosine[i] + 1) < teteps) cosine[i] = -1;
    }

    for (int i = 0; i < 6; i++) {
        angle[i] = REAL_ACOS(cosine[i]);
        sine[i]  = REAL_SIN(angle[i]);
        angle[i] /= twopi;
    }

    if constexpr (!compder) return;

    // Now compute derivatives of the angles with respect to the edge lengths
    // Since (see above):
    //                      num(i,j)
    // cos(ang(i,j)) = --------------------
    //                 sqrt(M(i,i)*M(j,j))

    //  d(ang(i,j))                      dnum(i,j)                             (M(i,i)dM(j,j) +M(j,j)*dM(i,i))
    // ------------sin(ang(i,j)) =  -----------------------   -0.5*num(i,j) -----------------------------------------
    //    dr(a,b)                   sqrt(M(i,i)M(j,j))dr(a,b)                  M(i,i)M(j,j) sqrt(M(i,i)M(j,j))

    // which we can rewrite as:

    // d(ang(i,j))                 cosine(i,j)    dnum(i,j)                    dM(j,j) +  dM(i,i))
    // ------------sin(ang(i,j)) = -----------  -----------  -0.5*cosine(i,j)( -------- + ---------)
    //     dr(a,b)                   num(i,j)       dr(a,b)                      M(j,j)      M(i,i)

    det[5] = det12; det[4] = det13; det[3] = det14;
    det[2] = det23; det[1] = det24; det[0] = det34;
    dist[0] = REAL_SQRT(r12sq); dist[1] = REAL_SQRT(r13sq); dist[2] = REAL_SQRT(r23sq);

    dminori[0][2] = -(val234 + 2*r24sq);
    dminori[1][1] = -(val134 + 2*r14sq);
    dminori[2][0] = -(val124 + 2*r14sq);
    dminori[3][0] = -(val123 + 2*r13sq); dminori[3][1] = -(val123 + 2*r12sq); dminori[3][2] = val123;

    dnum[5][0] = -2*val134+val123+val124; dnum[5][1] = 2*r12sq + val124; dnum[5][2] = -val124;
    dnum[4][0] = 2*r13sq+val134; dnum[4][1] = -2*val124 + val123 + val134; dnum[4][2] = -val134;
    dnum[3][0] = 2*r14sq+val134; dnum[3][1] = 2*r14sq + val124; dnum[3][2] = -2*r14sq;
    dnum[2][0] = 2*r23sq+val234; dnum[2][1] = -val234; dnum[2][2] = -2*val214 + val213 + val234;
    dnum[1][0] = 2*r24sq+val234; dnum[1][1] = -2*r24sq; dnum[1][2] = 2*r24sq + val214;
    dnum[0][0] = -2*r34sq; dnum[0][1] = 2*r34sq+val324; dnum[0][2] = 2*r34sq + val314;

    int k = 0;
    int jj;
    for (int i = 0; i < 3; i++) {
        for (int j = i+1; j < 4; j++) {
            jj = 5-k;
            if (det[k] != 0) {
                vala = cosine[jj]/sine[jj];
                val1 = -vala/det[k];
                val2 = vala/minori[j];
                val3 = vala/minori[i];
                for (int l = 0; l < 3; l++) {
                    deriv[jj][l] = val1*dnum[k][l]
                    +val2*dminori[j][l]+val3*dminori[i][l];
                    deriv[jj][l] *= 2*dist[l];
                }
            }
            else {
                vala = -val[i]*val[j]/sine[jj];
                for (int l = 0; l < 3; l++) {
                    deriv[jj][l] = vala*dnum[k][l];
                    deriv[jj][l] *= 2*dist[l];
                }
            }
            k++;
        }
    }
}

///////////////////////////////////////////////////////
//                                                   //
//  tet3dihedcos  --  compute three of siz dihedral  //
//                                                   //
///////////////////////////////////////////////////////

// "tet3dihedcos" computes three of the six dihedral angles
// of a tetrahedron from its edge lengths

template <bool compder>
inline void tet3dihedcos(real r12sq, real r13sq, real r14sq,
    real r23sq, real r24sq,real r34sq, real* cosine, real deriv[3][3])
{
    real val1, val2, val3, val4;
    real val123, val124, val134, val234;
    real val213, val214;
    real det12, det13, det23;
    real minori[4];
    real dminori[4][3] = {0};
    real dnum[3][3];
    real dist[3];

    // Define the Cayley Menger matrix:
    // M = ( 0      r12^2  r13^2  r14^2  1)
    //     ( r12^2  0      r23^2  r24^2  1)
    //     ( r13^2  r23^2  0      r34^2  1)
    //     ( r14^2  r24^2  r34^2  0      1)
    //     ( 1      1      1      1      0)
    // Compute all minors M(i,i): determinant of the Cayley-Menger matrix with row i
    // and column j removed
    // These determinants are of the form:
    // det = | 0 a b 1 |
    //       | a 0 c 1 |
    //       | b c 0 1 |
    //       | 1 1 1 0 |
    // then:
    // det = (c - a - b )^2 - 4ab

    val234 = (r34sq - r23sq - r24sq);
    val134 = (r34sq - r14sq - r13sq);
    val124 = (r24sq - r12sq - r14sq);
    val123 = (r23sq - r12sq - r13sq);

    minori[0] = val234*val234 - 4*r23sq*r24sq;
    minori[1] = val134*val134 - 4*r13sq*r14sq;
    minori[2] = val124*val124 - 4*r12sq*r14sq;
    minori[3] = val123*val123 - 4*r12sq*r13sq;

    val4 = 1.0/REAL_SQRT(-minori[0]);
    val3 = 1.0/REAL_SQRT(-minori[1]);
    val2 = 1.0/REAL_SQRT(-minori[2]);
    val1 = 1.0/REAL_SQRT(-minori[3]);

    // Now compute all angles (in fact, cosine of the angle):
    //           (-1)^(i+j) * det(Mij) 
    // cos(i,j)= ---------------------
    //            sqrt(M(i,i)*M(j,j))
    // where det(Mij) = M(i,j) is the determinant of the Cayley-Menger matrix with row i
    // and column j removed

    det12 = -2*r12sq*val134 - val123*val124;
    det13 = -2*r13sq*val124 - val123*val134;

    val213 = r13sq -r12sq -r23sq;
    val214 = r14sq -r12sq -r24sq;

    det23 = -2*r23sq*val214 - val213*val234;

    cosine[0] = det12*val1*val2;
    cosine[1] = det13*val1*val3;
    cosine[2] = det23*val1*val4;

    if constexpr (!compder) return;

    dminori[0][2] = -(val234 + 2*r24sq);
    dminori[1][1] = -(val134 + 2*r14sq);
    dminori[2][0] = -(val124 + 2*r14sq);
    dminori[3][0] = -(val123 + 2*r13sq);
    dminori[3][1] = -(val123 + 2*r12sq);
    dminori[3][2] = val123;

    dnum[0][0] = -2*val134+val123+val124;
    dnum[0][1] = 2*r12sq + val124;
    dnum[0][2] = -val124;

    dnum[1][0] = 2*r13sq + val134;
    dnum[1][1] = -2*val124 + val123 + val134;
    dnum[1][2] = -val134;

    dnum[2][0] = 2*r23sq + val234;
    dnum[2][1] = -val234;
    dnum[2][2] = -2*val214 + val213 + val234;

    dist[0] = REAL_SQRT(r12sq); dist[1] = REAL_SQRT(r13sq); dist[2] = REAL_SQRT(r23sq);

    for (int i = 0; i < 3; i++) {
        deriv[0][i] = dnum[0][i]*val1*val2 - cosine[0]*
        (dminori[2][i]/minori[2] + dminori[3][i]/minori[3]);
        deriv[1][i] = dnum[1][i]*val1*val3 - cosine[1]*
        (dminori[1][i]/minori[1] + dminori[3][i]/minori[3]);
        deriv[2][i] = dnum[2][i]*val1*val4 - cosine[2]*
        (dminori[0][i]/minori[0] + dminori[3][i]/minori[3]);
        deriv[0][i] *= 2*dist[i];
        deriv[1][i] *= 2*dist[i];
        deriv[2][i] *= 2*dist[i];
    }
}

/////////////////////////////////////////////////////////
//                                                     //
//  tetvorder  --  compute intersection of four balls  //
//                                                     //
/////////////////////////////////////////////////////////

// "tetvorder" computes the volume of the intersection of
// the tetrahedron formed by the center of 4 balls with
// the Voronoi cells corresponding to these balls

template <bool compder>
inline void tetvorder(real ra2,real rb2,real rc2,real rd2,
    real rab, real rac, real rad, real rbc, real rbd,
    real rcd, real rab2, real rac2, real rad2,real rbc2,
    real rbd2, real rcd2, real* cos_ang, real* sin_ang,
    real deriv[6][6], real& vola, real& volb, real& volc, 
    real& vold, real* dvola, real* dvolb, real* dvolc, real* dvold)
{
    real l1,l2,l3,l4,l5,l6;
    real val1,val2,val3,val4,val5,val6;
    real val1b,val2b,val3b,val4b,val5b,val6b;
    real val_ab,val_ac,val_bc,val_ad,val_bd,val_cd;
    real val1_ab,val1_ac,val1_ad,val1_bc,val1_bd,val1_cd;
    real val2_ab,val2_ac,val2_ad,val2_bc,val2_bd,val2_cd;
    real cos_abc,cos_acb,cos_bca,cos_abd,cos_adb,cos_bda;
    real cos_acd,cos_adc,cos_cda,cos_bcd,cos_bdc,cos_cdb;
    real rho_ab2,rho_ac2,rho_ad2,rho_bc2,rho_bd2,rho_cd2;
    real drho_ab2,drho_ac2,drho_ad2,drho_bc2,drho_bd2,drho_cd2;
    real dval1,dval2,dval3,dval4,dval5,dval6;
    real dval1b,dval2b,dval3b,dval4b,dval5b,dval6b;
    real cap_ab,cap_ac,cap_ad,cap_bc,cap_bd,cap_cd;
    real cosine_abc[3],cosine_abd[3],cosine_acd[3],cosine_bcd[3];
    real invsin[6],cotan[6];
    real deriv_abc[3][3],deriv_abd[3][3];
    real deriv_acd[3][3],deriv_bcd[3][3];
    real dinvsin[6][6],dcotan[6][6];
    real dval1_ab[6],dval1_ac[6],dval1_ad[6],dval1_bc[6];
    real dval1_bd[6],dval1_cd[6];
    real dval2_ab[6],dval2_ac[6],dval2_ad[6],dval2_bc[6];
    real dval2_bd[6],dval2_cd[6];
    real dcap_ab[6],dcap_ac[6],dcap_ad[6],dcap_bc[6];
    real dcap_bd[6],dcap_cd[6];

    l1 = plane_dist(ra2, rb2, rab2);
    l2 = plane_dist(ra2, rc2, rac2);
    l3 = plane_dist(ra2, rd2, rad2);
    l4 = plane_dist(rb2, rc2, rbc2);
    l5 = plane_dist(rb2, rd2, rbd2);
    l6 = plane_dist(rc2, rd2, rcd2);

    val1 = l1*rab; val2 = l2*rac; val3 = l3*rad;
    val4 = l4*rbc; val5 = l5*rbd; val6 = l6*rcd;

    val1b = rab-val1; val2b = rac-val2; val3b = rad-val3;
    val4b = rbc-val4; val5b = rbd-val5; val6b = rcd-val6;

    // We consider the tetrahedron (A,B,C,P_ABC) where P_ABC is the
    // point of intersection of the three spheres such that (A,B,C,P_ABC) is ccw.
    // The edge lengths in this tetrahedron are: rab, rac, rAP=ra, rbc, rBP=rb, rCP=rc

    tet3dihedcos<compder>(rab2, rac2, ra2, rbc2, rb2, rc2, cosine_abc, deriv_abc);

    // repeat for tetrahedron (A,B,D,P_ABD)
    tet3dihedcos<compder>(rab2, rad2, ra2, rbd2, rb2, rd2, cosine_abd, deriv_abd);

    // repeat for tetrahedron (A,C,D,P_ACD)
    tet3dihedcos<compder>(rac2, rad2, ra2, rcd2, rc2, rd2, cosine_acd, deriv_acd);

    // repeat for tetrahedron (B,C,D,P_BCD)
    tet3dihedcos<compder>(rbc2, rbd2, rb2, rcd2, rc2, rd2, cosine_bcd, deriv_bcd);

    cos_abc = cosine_abc[0]; cos_acb = cosine_abc[1]; cos_bca = cosine_abc[2];
    cos_abd = cosine_abd[0]; cos_adb = cosine_abd[1]; cos_bda = cosine_abd[2];
    cos_acd = cosine_acd[0]; cos_adc = cosine_acd[1]; cos_cda = cosine_acd[2];
    cos_bcd = cosine_bcd[0]; cos_bdc = cosine_bcd[1]; cos_cdb = cosine_bcd[2];

    rho_ab2 = ra2 - val1b*val1b; rho_ac2 = ra2 - val2b*val2b;
    rho_ad2 = ra2 - val3b*val3b; rho_bc2 = rb2 - val4b*val4b;
    rho_bd2 = rb2 - val5b*val5b; rho_cd2 = rc2 - val6b*val6b;

    for (int i = 0; i < 6; i++) {
        if (REAL_ABS(sin_ang[i]) < teteps) {
            invsin[i] = 0;
            cotan[i] = 0;
        }
        else {
            invsin[i] = 1.0/sin_ang[i];
            cotan[i] = cos_ang[i]*invsin[i];
        }
    }

    val_ab = -(cos_abc*cos_abc+cos_abd*cos_abd)*cotan[0]
        + 2*cos_abc*cos_abd*invsin[0];
    val_ac = -(cos_acb*cos_acb+cos_acd*cos_acd)*cotan[1]
        + 2*cos_acb*cos_acd*invsin[1];
    val_ad = -(cos_adb*cos_adb+cos_adc*cos_adc)*cotan[2]
        + 2*cos_adb*cos_adc*invsin[2];
    val_bc = -(cos_bca*cos_bca+cos_bcd*cos_bcd)*cotan[3]
        + 2*cos_bca*cos_bcd*invsin[3];
    val_bd = -(cos_bda*cos_bda+cos_bdc*cos_bdc)*cotan[4]
        + 2*cos_bda*cos_bdc*invsin[4];
    val_cd = -(cos_cda*cos_cda+cos_cdb*cos_cdb)*cotan[5]
        + 2*cos_cda*cos_cdb*invsin[5];

    cap_ab = rho_ab2*val_ab; cap_ac = rho_ac2*val_ac;
    cap_ad = rho_ad2*val_ad; cap_bc = rho_bc2*val_bc;
    cap_bd = rho_bd2*val_bd; cap_cd = rho_cd2*val_cd;

    vola = (val1b*cap_ab+val2b*cap_ac+val3b*cap_ad)/6;
    volb = (val1*cap_ab+val4b*cap_bc+val5b*cap_bd)/6;
    volc = (val2*cap_ac+val4*cap_bc+val6b*cap_cd)/6;
    vold = (val3*cap_ad+val5*cap_bd+val6*cap_cd)/6;

    if constexpr (!compder) return;

    dval1b = l1; dval2b = l2; dval3b = l3;
    dval4b = l4; dval5b = l5; dval6b = l6;

    dval1 = 1-l1; dval2 = 1-l2; dval3 = 1-l3;
    dval4 = 1-l4; dval5 = 1-l5; dval6 = 1-l6;

    drho_ab2 = -2*dval1b*val1b; drho_ac2 = -2*dval2b*val2b;
    drho_ad2 = -2*dval3b*val3b; drho_bc2 = -2*dval4b*val4b;
    drho_bd2 = -2*dval5b*val5b; drho_cd2 = -2*dval6b*val6b;

    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            dcotan[i][j] = -deriv[i][j]*(1.+cotan[i]*cotan[i]);
            dinvsin[i][j] = -deriv[i][j]*cotan[i]*invsin[i];
        }
    }

    val1_ab = cos_abc*cos_abc+cos_abd*cos_abd;
    val2_ab = 2*cos_abc*cos_abd;

    dval1_ab[0] = 2*(deriv_abc[0][0]*cos_abc+deriv_abd[0][0]*cos_abd);
    dval1_ab[1] = 2*deriv_abc[0][1]*cos_abc;
    dval1_ab[2] = 2*deriv_abd[0][1]*cos_abd;
    dval1_ab[3] = 2*deriv_abc[0][2]*cos_abc;
    dval1_ab[4] = 2*deriv_abd[0][2]*cos_abd;
    dval1_ab[5] = 0;

    dval2_ab[0] = 2*(deriv_abc[0][0]*cos_abd+deriv_abd[0][0]*cos_abc);
    dval2_ab[1] = 2*deriv_abc[0][1]*cos_abd;
    dval2_ab[2] = 2*deriv_abd[0][1]*cos_abc;
    dval2_ab[3] = 2*deriv_abc[0][2]*cos_abd;
    dval2_ab[4] = 2*deriv_abd[0][2]*cos_abc;
    dval2_ab[5] = 0;

    for (int i = 0; i < 6; i++) {
        dcap_ab[i] = - dval1_ab[i]*cotan[0] - val1_ab*dcotan[0][i]
            + dval2_ab[i]*invsin[0] + val2_ab*dinvsin[0][i];
        dcap_ab[i] = rho_ab2*dcap_ab[i];
    }
    dcap_ab[0] += drho_ab2*val_ab;

    val1_ac = cos_acb*cos_acb + cos_acd*cos_acd;
    val2_ac = 2*cos_acb*cos_acd;

    dval1_ac[0] = 2*deriv_abc[1][0]*cos_acb;
    dval1_ac[1] = 2*(deriv_abc[1][1]*cos_acb+deriv_acd[0][0]*cos_acd);
    dval1_ac[2] = 2*deriv_acd[0][1]*cos_acd;
    dval1_ac[3] = 2*deriv_abc[1][2]*cos_acb;
    dval1_ac[4] = 0;
    dval1_ac[5] = 2*deriv_acd[0][2]*cos_acd;

    dval2_ac[0] = 2*deriv_abc[1][0]*cos_acd;
    dval2_ac[1] = 2*(deriv_abc[1][1]*cos_acd+deriv_acd[0][0]*cos_acb);
    dval2_ac[2] = 2*deriv_acd[0][1]*cos_acb;
    dval2_ac[3] = 2*deriv_abc[1][2]*cos_acd;
    dval2_ac[4] = 0;
    dval2_ac[5] = 2*deriv_acd[0][2]*cos_acb;

    for (int i = 0; i < 6; i++) {
        dcap_ac[i] = -dval1_ac[i]*cotan[1] - val1_ac*dcotan[1][i]
            + dval2_ac[i]*invsin[1] + val2_ac*dinvsin[1][i];
        dcap_ac[i] = rho_ac2*dcap_ac[i];
    }
    dcap_ac[1] += drho_ac2*val_ac;

    val1_ad = cos_adb*cos_adb + cos_adc*cos_adc;
    val2_ad = 2*cos_adb*cos_adc;

    dval1_ad[0] = 2*deriv_abd[1][0]*cos_adb;
    dval1_ad[1] = 2*deriv_acd[1][0]*cos_adc;
    dval1_ad[2] = 2*(deriv_abd[1][1]*cos_adb+deriv_acd[1][1]*cos_adc);
    dval1_ad[3] = 0;
    dval1_ad[4] = 2*deriv_abd[1][2]*cos_adb;
    dval1_ad[5] = 2*deriv_acd[1][2]*cos_adc;

    dval2_ad[0] = 2*deriv_abd[1][0]*cos_adc;
    dval2_ad[1] = 2*deriv_acd[1][0]*cos_adb;
    dval2_ad[2] = 2*(deriv_abd[1][1]*cos_adc+deriv_acd[1][1]*cos_adb);
    dval2_ad[3] = 0;
    dval2_ad[4] = 2*deriv_abd[1][2]*cos_adc;
    dval2_ad[5] = 2*deriv_acd[1][2]*cos_adb;

    for (int i = 0; i < 6; i++) {
        dcap_ad[i] = -dval1_ad[i]*cotan[2] - val1_ad*dcotan[2][i]
            + dval2_ad[i]*invsin[2] + val2_ad*dinvsin[2][i];
        dcap_ad[i] = rho_ad2*dcap_ad[i];
    }
    dcap_ad[2] += drho_ad2*val_ad;

    val1_bc = cos_bca*cos_bca + cos_bcd*cos_bcd;
    val2_bc = 2*cos_bca*cos_bcd;

    dval1_bc[0] = 2*deriv_abc[2][0]*cos_bca;
    dval1_bc[1] = 2*deriv_abc[2][1]*cos_bca;
    dval1_bc[2] = 0;
    dval1_bc[3] = 2*(deriv_abc[2][2]*cos_bca+deriv_bcd[0][0]*cos_bcd);
    dval1_bc[4] = 2*deriv_bcd[0][1]*cos_bcd;
    dval1_bc[5] = 2*deriv_bcd[0][2]*cos_bcd;

    dval2_bc[0] = 2*deriv_abc[2][0]*cos_bcd;
    dval2_bc[1] = 2*deriv_abc[2][1]*cos_bcd;
    dval2_bc[2] = 0;
    dval2_bc[3] = 2*(deriv_abc[2][2]*cos_bcd+deriv_bcd[0][0]*cos_bca);
    dval2_bc[4] = 2*deriv_bcd[0][1]*cos_bca;
    dval2_bc[5] = 2*deriv_bcd[0][2]*cos_bca;

    for (int i = 0; i < 6; i++) {
        dcap_bc[i] = -dval1_bc[i]*cotan[3] - val1_bc*dcotan[3][i]
            + dval2_bc[i]*invsin[3] + val2_bc*dinvsin[3][i];
        dcap_bc[i] = rho_bc2*dcap_bc[i];
    }
    dcap_bc[3] += drho_bc2*val_bc;

    val1_bd = cos_bda*cos_bda + cos_bdc*cos_bdc;
    val2_bd = 2*cos_bda*cos_bdc;

    dval1_bd[0] = 2*deriv_abd[2][0]*cos_bda;
    dval1_bd[1] = 0;
    dval1_bd[2] = 2*deriv_abd[2][1]*cos_bda;
    dval1_bd[3] = 2*deriv_bcd[1][0]*cos_bdc;
    dval1_bd[4] = 2*(deriv_abd[2][2]*cos_bda+deriv_bcd[1][1]*cos_bdc);
    dval1_bd[5] = 2*deriv_bcd[1][2]*cos_bdc;

    dval2_bd[0] = 2*deriv_abd[2][0]*cos_bdc;
    dval2_bd[1] = 0;
    dval2_bd[2] = 2*deriv_abd[2][1]*cos_bdc;
    dval2_bd[3] = 2*deriv_bcd[1][0]*cos_bda;
    dval2_bd[4] = 2*(deriv_abd[2][2]*cos_bdc+deriv_bcd[1][1]*cos_bda);
    dval2_bd[5] = 2*deriv_bcd[1][2]*cos_bda;

    for (int i = 0; i < 6; i++) {
        dcap_bd[i] = -dval1_bd[i]*cotan[4] - val1_bd*dcotan[4][i]
            + dval2_bd[i]*invsin[4] + val2_bd*dinvsin[4][i];
        dcap_bd[i] = rho_bd2*dcap_bd[i];
    }
    dcap_bd[4] += drho_bd2*val_bd;

    val1_cd = cos_cda*cos_cda + cos_cdb*cos_cdb;
    val2_cd = 2*cos_cda*cos_cdb;

    dval1_cd[0] = 0;
    dval1_cd[1] = 2*deriv_acd[2][0]*cos_cda;
    dval1_cd[2] = 2*deriv_acd[2][1]*cos_cda;
    dval1_cd[3] = 2*deriv_bcd[2][0]*cos_cdb;
    dval1_cd[4] = 2*deriv_bcd[2][1]*cos_cdb;
    dval1_cd[5] = 2*(deriv_acd[2][2]*cos_cda+deriv_bcd[2][2]*cos_cdb);

    dval2_cd[0] = 0;
    dval2_cd[1] = 2*deriv_acd[2][0]*cos_cdb;
    dval2_cd[2] = 2*deriv_acd[2][1]*cos_cdb;
    dval2_cd[3] = 2*deriv_bcd[2][0]*cos_cda;
    dval2_cd[4] = 2*deriv_bcd[2][1]*cos_cda;
    dval2_cd[5] = 2*(deriv_acd[2][2]*cos_cdb+deriv_bcd[2][2]*cos_cda);

    for (int i = 0; i < 6; i++) {
        dcap_cd[i] = -dval1_cd[i]*cotan[5]-val1_cd*dcotan[5][i]
            +dval2_cd[i]*invsin[5]+val2_cd*dinvsin[5][i];
        dcap_cd[i] = rho_cd2*dcap_cd[i];
    }
    dcap_cd[5] += drho_cd2*val_cd;

    for (int i = 0; i < 6; i++) {
        dvola[i] = (val1b*dcap_ab[i]+val2b*dcap_ac[i]
            + val3b*dcap_ad[i])/6;
        dvolb[i] = (val1*dcap_ab[i]+val4b*dcap_bc[i]
            + val5b*dcap_bd[i])/6;
        dvolc[i] = (val2*dcap_ac[i]+val4*dcap_bc[i]
            + val6b*dcap_cd[i])/6;
        dvold[i] = (val3*dcap_ad[i]+val5*dcap_bd[i]
            + val6*dcap_cd[i])/6;
    }

    dvola[0] += dval1b*cap_ab/6;
    dvola[1] += dval2b*cap_ac/6;
    dvola[2] += dval3b*cap_ad/6;
    dvolb[0] += dval1*cap_ab/6;
    dvolb[3] += dval4b*cap_bc/6;
    dvolb[4] += dval5b*cap_bd/6;
    dvolc[1] += dval2*cap_ac/6;
    dvolc[3] += dval4*cap_bc/6;
    dvolc[5] += dval6b*cap_cd/6;
    dvold[2] += dval3*cap_ad/6;
    dvold[4] += dval5*cap_bd/6;
    dvold[5] += dval6*cap_cd/6;
}
}
