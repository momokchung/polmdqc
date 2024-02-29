// Author: Moses KJ Chung
// Year:   2024

#pragma once
#include "alfc.h"
#include "libfunc.h"

namespace polmdqc
{
////////////////////////////////////////////////////////////////
//                                                            //
//  alftrig  --  checks if triangle belongs to alpha complex  //
//                                                            //
////////////////////////////////////////////////////////////////

// "alftrig" checks if a triangle ABC belons to the alpha complex

// triattach checks if triangles are attached
inline void triattach(real* a, real* b, real* c, real* d,
    real ra, real rb, real rc, real rd, real S[3][4], real T[2][3],
    real Dabc, int& testa, int& memory)
{
    testa = 0;

    // We need to compute:
    // Det1 = Minor(a,b,c,d,2,3,4,0)
    // Det2 = Minor(a,b,c,d,1,3,4,0)
    // Det3 = Minor(a,b,c,d,1,2,4,0)
    // Deter= Minor(a,b,c,d,1,2,3,0)

    real Det1 = -d[1]*S[2][3] + d[2]*S[1][3] - d[3]*S[1][2] + T[1][2];
    real Det2 = -d[0]*S[2][3] + d[2]*S[0][3] - d[3]*S[0][2] + T[0][2];
    real Det3 = -d[0]*S[1][3] + d[1]*S[0][3] - d[3]*S[0][1] + T[0][1];
    real Deter = -d[0]*S[1][2] + d[1]*S[0][2] - d[2]*S[0][1] + Dabc;

    // Check if the face is "attached" to the fourth
    // vertex of the parent tetrahedron

    real test = Det1*S[1][2]+Det2*S[0][2]+Det3*S[0][1]-2*Deter*Dabc;

    // check for problems, in which case change precision

    // if (REAL_ABS(test) < alfeps) {
    //     test = 0.;
    //     memory = 1;
    // }

    // if no problem, set testa to true if test > 0
    if (test > 0) testa = 1;
    return;
}

// triradius computes the radius of the
// smallest circumsphere to a triangle
inline void triradius(real* a, real* b, real* c, real ra, real rb,
    real rc, real S[3][4], real T[2][3], real Dabc, int& testr, real alpha, int& memory)
{
    testr = 0;
    real sums2 = S[0][1]*S[0][1] + S[0][2]*S[0][2] + S[1][2]*S[1][2];
    real d0 = sums2;
    real d1 = S[0][2]*S[2][3] + S[0][1]*S[1][3] - 2*Dabc*S[1][2];
    real d2 = S[0][1]*S[0][3] - S[1][2]*S[2][3] - 2*Dabc*S[0][2];
    real d3 = S[1][2]*S[1][3] + S[0][2]*S[0][3] + 2*Dabc*S[0][1];
    real d4 = S[0][1]*T[0][1] + S[0][2]*T[0][2] + S[1][2]*T[1][2] - 2*Dabc*Dabc;
    real num  = 4*(d1*d1+d2*d2+d3*d3) + 16*d0*d4;

    // if (REAL_ABS(alpha-num) < alfeps) {
    //     num = alpha;
    //     memory = 1;
    // }

    if (alpha > num) testr = 1;
    return;
}

inline void alftrig(real* a, real* b, real* c, real* d, real* e,
    real ra,real rb, real rc, real rd, real re, int ie, 
    int& irad,int& iattach, real alpha)
{
    real Dabc;
    real Sab[3][4],Sac[3][4],Sbc[3][4];
    real S[3][4],T[2][3];

    iattach = 0;
    irad = 0;

    real val = a[3]+b[3] -2*(a[0]*b[0]+a[1]*b[1]+a[2]*b[2]+ra*rb);
    if (val > 0) return;
    val = a[3]+c[3] -2*(a[0]*c[0]+a[1]*c[1]+a[2]*c[2]+ra*rc);
    if (val > 0) return;
    val = b[3]+c[3] -2*(b[0]*c[0]+b[1]*c[1]+b[2]*c[2]+rb*rc);
    if (val > 0) return;

    // Perform computation in floating points; if a problem occurs,
    // switch precision

    // 1. Computes all Minors Smn(i,j) = M(m,n,i,j) = Det | m(i)  m(j) |
    //                                                    | n(i)  n(j) |
    // m,n are two vertices of the triangle, i and j correspond
    // to two of the coordinates of the vertices

    // for all i in [0,2] and all j in [i+1,3]
    for (int i = 0; i < 3; i++) {
        for (int j = i+1; j < 4; j++) {
            Sab[i][j] = a[i]*b[j]-a[j]*b[i];
            Sac[i][j] = a[i]*c[j]-a[j]*c[i];
            Sbc[i][j] = b[i]*c[j]-b[j]*c[i];
        }
    }

    // Now compute all Minors 
    //     S(i,j) = M(a,b,c,i,j,0) = Det | a(i) a(j) 1 |
    //                                   | b(i) b(j) 1 |
    //                                   | c(i) c(j) 1 |
    // a,b,c are the 3 vertices of the triangle, i and j correspond
    // to two of the coordinates of the vertices

    // for all i in [0,2] and all j in [i+1,3]
    for (int i = 0; i < 3; i++) {
        for (int j = i+1; j < 4; j++) {
            S[i][j] = Sbc[i][j] - Sac[i][j] + Sab[i][j];
        }
    }

    // Now compute all Minors
    //     T(i,j) = M(a,b,c,i,j,4) = Det | a(i) a(j) a(4) |
    //                                   | b(i) b(j) b(4) |
    //                                   | c(i) c(j) c(4) |

    // for all i in [0,1] and all j in [i+1,2]
    for (int i = 0; i < 2; i++) {
        for (int j = i+1; j < 3; j++) {
            T[i][j] = a[3]*Sbc[i][j] - b[3]*Sac[i][j] + c[3]*Sab[i][j];
        }
    }

    // Finally,  need Dabc = M(a,b,c,1,2,3) = Det | a(1) a(2) a(3) |
    //                                            | b(2) b(2) b(3) |
    //                                            | c(3) c(2) c(3) |

    Dabc = a[0]*Sbc[1][2] - b[0]*Sac[1][2] + c[0]*Sab[1][2];

    // first check if a,b,c attached to d:
    int memory = 0;
    int attach;
    triattach(a, b, c, d, ra, rb, rc, rd, S, T, Dabc, attach, memory);

    // If attached, we can stop there, the triangle will 
    // not be part of the alpha complex
    if (attach == 1) {
        iattach = 1;
        return;
    }

    // if e exists, check if a,b,c attached to e:
    if (ie >= 0) {
        triattach(a, b, c, e, ra, rb, rc, re, S, T, Dabc, attach, memory);

        // If attached, we can stop there, the triangle will
        // not be part of the alpha complex
        if (attach == 1) {
            iattach = 1;
            return;
        }
    }

    // Now check if alpha is bigger than the radius of the sphere orthogonal
    // to the three balls at A, B, C:

    int testr;
    triradius(a, b, c, ra, rb, rc, S, T, Dabc, testr, alpha, memory);

    if (testr == 1) irad = 1;
}
}
