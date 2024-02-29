// Author: Moses KJ Chung
// Year:   2024

#pragma once
#include "alfc.h"
#include "alffunc.h"
#include "alftetra.h"
#include "libfunc.h"
#include "precision.h"

namespace polmdqc
{
//////////////////////////////////////////////////////////
//                                                      //
//  alftetra  --  computes radius of orthogonal sphere  //
//                                                      //
//////////////////////////////////////////////////////////

// "alftetra" computes the radius R of the sphere orthogonal
// to the four spheres that define a tetrahedron [A,B,C,D]

// // computes the radius R of the circumsphere containing a tetrahedron [A,B,C,D]
// inline void tetrad(real* a, real* b, real* c, real* d, real ra,
//     real rb, real rc, real rd, int& testr, real alpha)
// {
//     int k,coef;
//     real wa,wb,wc,wd;
//     real temp1,temp2,temp3;
//     real det1,det2,det3,det4;
//     real Dabc,Dabd,Dacd,Dbcd,Dabcd;
//     real den,num;
//     real Sab[3],Sac[3],Sad[3],Sbc[3],Sbd[3],Scd[3];
//     real Sa[3],Sb[3],Sc[3],Sd[3];
//     real Sam1[3],Sbm1[3],Scm1[3],Sdm1[3];
//     real Deter[3];

//     buildweight(a[0], a[1], a[2], ra, wa);
//     buildweight(b[0], b[1], b[2], rb, wb);
//     buildweight(c[0], c[1], c[2], rc, wc);
//     buildweight(d[0], d[1], d[2], rd, wd);

//     // 1. Computes all Minors Smn(i+j-2)= M(m,n,i,j) = Det | m(i)  m(j) |
//     //                                                     | n(i)  n(j) |
//     // for all i in [1,2] and all j in [i+1,3]

//     for (int i = 0; i < 2; i++) {
//         for (int j = i+1; j < 3 ; j++) {
//             k = i + j - 1;
//             temp1 = a[j] * b[i];
//             temp2 = a[i] * b[j];
//             Sab[k] = psub(temp2,temp1);
//             temp1 = a[j] * c[i];
//             temp2 = a[i] * c[j];
//             Sac[k] = psub(temp2,temp1);
//             temp1 = a[j] * d[i];
//             temp2 = a[i] * d[j];
//             Sad[k] = psub(temp2,temp1);
//             temp1 = b[j] * c[i];
//             temp2 = b[i] * c[j];
//             Sbc[k] = psub(temp2,temp1);
//             temp1 = b[j] * d[i];
//             temp2 = b[i] * d[j];
//             Sbd[k] = psub(temp2,temp1);
//             temp1 = c[j] * d[i];
//             temp2 = c[i] * d[j];
//             Scd[k] = psub(temp2,temp1);
//         }
//     }

//     // Now compute all Minors 
//     //     Sq(i+j-2) = M(m,n,p,i,j,0) = Det | m(i) m(j) 1 |
//     //                                      | n(i) n(j) 1 |
//     //                                      | p(i) p(j) 1 |

//     // and all Minors
//     //     Det(i+j-2) = M(m,n,p,q,i,j,4,0) = Det | m(i) m(j) m(4) 1 |
//     //                                           | n(i) n(j) n(4) 1 |
//     //                                           | p(i) p(j) p(4) 1 |
//     //                                           | q(i) q(j) q(4) 1 |

//     // m,n,p,q are the four vertices of the tetrahedron, i and j correspond
//     // to two of the coordinates of the vertices, and m(4) refers to the
//     // "weight" of vertices m

//     for (int i = 0; i < 3; i++) {
//         temp1 = psub(Scd[i],Sbd[i]);
//         Sa[i] = padd(temp1,Sbc[i]);
//         temp2 = Sa[i] * wa;
//         temp1 = psub(Scd[i],Sad[i]);
//         Sb[i] = padd(temp1,Sac[i]);
//         temp3 = Sb[i] * wb;
//         temp2 = psub(temp2,temp3);
//         temp1 = psub(Sbd[i],Sad[i]);
//         Sc[i] = padd(temp1,Sab[i]);
//         temp3 = Sc[i] * wc;
//         temp2 = padd(temp2,temp3);
//         temp1 = psub(Sbc[i],Sac[i]);
//         Sd[i] = padd(temp1,Sab[i]);
//         temp3 = Sd[i] * wd;
//         Deter[i] = psub(temp2,temp3);
//         Sam1[i] = -Sa[i];
//         Sbm1[i] = -Sb[i];
//         Scm1[i] = -Sc[i];
//         Sdm1[i] = -Sd[i];
//     }

//     // Now compute the determinant needed to compute the radius of the
//     // circumsphere of the tetrahedron :
//     //     Det1 = Minor(a,b,c,d,4,2,3,0)
//     //     Det2 = Minor(a,b,c,d,1,3,4,0)
//     //     Det3 = Minor(a,b,c,d,1,2,4,0)
//     //     Det4 = Minor(a,b,c,d,1,2,3,0)

//     det1 = Deter[2];
//     det2 = Deter[1];
//     det3 = Deter[0];

//     temp1 = a[0] * Sa[2];
//     temp2 = b[0] * Sb[2];
//     temp3 = psub(temp1,temp2);
//     temp1 = c[0] * Sc[2];
//     temp2 = d[0] * Sd[2];
//     temp1 = psub(temp1,temp2);
//     det4 = padd(temp1,temp3);

//     // Now compute all minors:
//     //     Dmnp = Minor(m, n, p, 1, 2, 3) = Det | m(1) m(2) m(3) |
//     //                                          | n(1) n(2) n(3) |
//     //                                          | p(1) p(2) p(3) |

//     temp1 = a[0] * Sbc[2];
//     temp2 = b[0] * Sac[2];
//     temp3 = psub(temp1,temp2);
//     temp1 = c[0] * Sab[2];
//     Dabc = padd(temp3,temp1);

//     temp1 = a[0] * Sbd[2];
//     temp2 = b[0] * Sad[2];
//     temp3 = psub(temp1,temp2);
//     temp1 = d[0] * Sab[2];
//     Dabd = padd(temp3,temp1);

//     temp1 = a[0] * Scd[2];
//     temp2 = c[0] * Sad[2];
//     temp3 = psub(temp1,temp2);
//     temp1 = d[0] * Sac[2];
//     Dacd = padd(temp3,temp1);

//     temp1 = b[0] * Scd[2];
//     temp2 = c[0] * Sbd[2];
//     temp3 = psub(temp1,temp2);
//     temp1 = d[0] * Sbc[2];
//     Dbcd = padd(temp3,temp1);

//     // We also need :
//     //     Det = Det | m(1) m(2) m(3) m(4) |
//     //               | n(1) n(2) n(3) n(4) |
//     //               | p(1) p(2) p(3) p(4) |
//     //               | q(1) q(2) q(3) q(4) |

//     temp1 = wa * Dbcd;
//     temp2 = wb * Dacd;
//     temp3 = psub(temp2,temp1);
//     temp1 = wc * Dabd;
//     temp2 = wd * Dabc;
//     temp1 = psub(temp2,temp1);
//     Dabcd = padd(temp3,temp1);

//     // The radius of the circumsphere of the weighted tetrahedron is then:
//     // r_t = (Det1*Det1 + Det2*Det2 + Det3*Det3 + 4*Det4*Dabcd)/(4*Det4*Det4)
//     coef=4;
//     temp1 = det4 * det4;
//     den = temp1 * coef;

//     temp1 = det1 * det1;
//     temp2 = det2 * det2;
//     temp1 = padd(temp1,temp2);
//     temp2 = det3 * det3;
//     temp1 = padd(temp1,temp2);
//     temp2 = det4 * Dabcd;
//     temp2 = temp2 * coef;
//     num = padd(temp1,temp2);

//     temp1 = den * alpha;
//     temp2 = psub(num,temp1);

//     // If tetrahedron is part of the alpha shape, then the 4 triangles,
//     // the 6 edges and the four vertices are also part of the alpha
//     // complex
//     int temp2sgn = sgn(temp2);
//     if(!(temp2sgn > 0)) 
//     {
//         testr = 1;
//     }
//     else
//     {
//         testr = 0;
//     }
// }

inline void alftetra(real* a, real* b, real* c, real* d,
    real ra, real rb, real rc, real rd, int& iflag, real alpha)
{
    real D1, D2, D3, D4, Det;
    real Dabc, Dabd, Dacd, Dbcd;
    real num,den;
    real test,val;
    real Sab[3], Sac[3], Sad[3], Sbc[3], Sbd[3], Scd[3];
    real Sa[3], Sb[3], Sc[3], Sd[3];
    real Deter[3];

    iflag = 0;
    val = a[3]+b[3] -2*(a[0]*b[0]+a[1]*b[1]+a[2]*b[2]+ra*rb);
    if (val > 0) return;
    val = a[3]+c[3] -2*(a[0]*c[0]+a[1]*c[1]+a[2]*c[2]+ra*rc);
    if (val > 0) return;
    val = a[3]+d[3] -2*(a[0]*d[0]+a[1]*d[1]+a[2]*d[2]+ra*rd);
    if (val > 0) return;
    val = b[3]+c[3] -2*(b[0]*c[0]+b[1]*c[1]+b[2]*c[2]+rb*rc);
    if (val > 0) return;
    val = b[3]+d[3] -2*(b[0]*d[0]+b[1]*d[1]+b[2]*d[2]+rb*rd);
    if (val > 0) return;
    val = c[3]+d[3] -2*(c[0]*d[0]+c[1]*d[1]+c[2]*d[2]+rc*rd);
    if (val > 0) return;

    // Perform computation in floating points; if a problem occurs, switch precision

    // 1. Computes all Minors Smn(i+j-2)= M(m,n,i,j) = Det | m(i)  m(j) |
    //                                                     | n(i)  n(j) |
    //    for all i in [0,1] and all j in [i+1,2]

    for (int i = 0; i < 2; i++) {
        for (int j = i+1; j < 3; j++) {
            int k = i+j-1;
            Sab[k] = a[i]*b[j]-a[j]*b[i];
            Sac[k] = a[i]*c[j]-a[j]*c[i];
            Sad[k] = a[i]*d[j]-a[j]*d[i];
            Sbc[k] = b[i]*c[j]-b[j]*c[i];
            Sbd[k] = b[i]*d[j]-b[j]*d[i];
            Scd[k] = c[i]*d[j]-c[j]*d[i];
        }
    }

    // Now compute all Minors 
    //     Sq(i+j-2) = M(m,n,p,i,j,0) = Det | m(i) m(j) 1 |
    //                                      | n(i) n(j) 1 |
    //                                      | p(i) p(j) 1 |

    // and all Minors
    //     Det(i+j-2) = M(m,n,p,q,i,j,4,0) = Det | m(i) m(j) m(4) 1 |
    //                                           | n(i) n(j) n(4) 1 |
    //                                           | p(i) p(j) p(4) 1 |
    //                                           | q(i) q(j) q(4) 1 |

    // m,n,p,q are the four vertices of the tetrahedron, i and j correspond
    // to two of the coordinates of the vertices, and m(4) refers to the
    // "weight" of vertices m

    for (int i = 0; i < 3; i++) {
        Sa[i] = Scd[i] - Sbd[i] + Sbc[i];
        Sb[i] = Scd[i] - Sad[i] + Sac[i];
        Sc[i] = Sbd[i] - Sad[i] + Sab[i];
        Sd[i] = Sbc[i] - Sac[i] + Sab[i];
    }

    for (int i = 0; i < 3; i++) {
        Deter[i] = a[3]*Sa[i]-b[3]*Sb[i]+c[3]*Sc[i]-d[3]*Sd[i];
    }

    // Now compute the determinant needed to compute the radius of the
    // sphere orthogonal to the four balls that define the tetrahedron :
    //     D1 = Minor(a,b,c,d,4,2,3,0)
    //     D2 = Minor(a,b,c,d,1,3,4,0)
    //     D3 = Minor(a,b,c,d,1,2,4,0)
    //     D4 = Minor(a,b,c,d,1,2,3,0)

    D1 = Deter[2];
    D2 = Deter[1];
    D3 = Deter[0];
    D4 = a[0]*Sa[2]-b[0]*Sb[2]+c[0]*Sc[2]-d[0]*Sd[2];

    // Now compute all minors:
    //     Dmnp = Minor(m,n,p,1,2,3) = Det | m(1) m(2) m(3) |
    //                                     | n(1) n(2) n(3) |
    //                                     | p(1) p(2) p(3) |

    Dabc = a[0]*Sbc[2]-b[0]*Sac[2] + c[0]*Sab[2];
    Dabd = a[0]*Sbd[2]-b[0]*Sad[2] + d[0]*Sab[2];
    Dacd = a[0]*Scd[2]-c[0]*Sad[2] + d[0]*Sac[2];
    Dbcd = b[0]*Scd[2]-c[0]*Sbd[2] + d[0]*Sbc[2];

    // We also need :
    //     Det = Det | m(1) m(2) m(3) m(4) |
    //               | n(1) n(2) n(3) n(4) |
    //               | p(1) p(2) p(3) p(4) |
    //               | q(1) q(2) q(3) q(4) |

    Det = -a[3]*Dbcd + b[3]*Dacd -c[3]*Dabd + d[3]*Dabc;

    // the radius of the circumsphere of the weighted tetrahedron is then:

    num = D1*D1 + D2*D2 + D3*D3 + 4*D4*Det;
    den = 4*D4*D4;

    // if this radius is too close to the value of ALPHA, we switch precision
    test = alpha*den - num;
    // int itest;
    // if (REAL_ABS(test) < alfeps) {
    //     tetrad(a, b, c, d, ra, rb, rc, rd, itest, alpha);
    //     test = itest;
    // }

    // The spectrum for a tetrahedron is [R_t Infinity[. If ALPHA is in
    // that interval, the tetrahedron is part of the alpha shape, otherwise
    // it is discarded
    // If tetrahedron is part of the alpha shape, then the 4 triangles,
    // the 6 edges and the four vertices are also part of the alpha
    // complex

    iflag = 0;
    if (test > 0) iflag = 1;
}
}
