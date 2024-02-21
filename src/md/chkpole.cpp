// Author: Moses KJ Chung
// Year:   2023

#include "atoms.h"
#include "chkpole.h"
#include "mpole.h"
#include "repel.h"
#include <cmath>

namespace polmdqc
{
/////////////////////////////////////////////////////
//                                                 //
//  chkpole  --  check multipoles at chiral sites  //
//                                                 //
/////////////////////////////////////////////////////

// "chkpole" inverts multipole moments as necessary at atoms
// with chiral local reference frame definitions

inline void chkpoleAtomI(int i, const MDQCArray<real> x, const MDQCArray<real> y, const MDQCArray<real> z, const MDQCArray<LocalFrame> polaxe,
                         const MDQCArray<int> xaxis, MDQCArray<int> yaxis, const MDQCArray<int> zaxis, MDQCArray2D<real,maxpole> pole);

void chkpole()
{
    // loop over multipoles and test for chirality inversion
    #pragma omp parallel for default(private)   \
        shared(n,x,y,z,polaxe,xaxis,yaxis,zaxis,pole)
    for (int i = 0; i < n; i++) {
        chkpoleAtomI(i, x, y, z, polaxe, xaxis, yaxis, zaxis, pole);
    }
}

void chkrepole()
{
    // loop over multipoles and test for chirality inversion
    #pragma omp parallel for default(private)   \
        shared(n,x,y,z,polaxe,xaxis,yaxis,zaxis,repole)
    for (int i = 0; i < n; i++) {
        chkpoleAtomI(i, x, y, z, polaxe, xaxis, yaxis, zaxis, repole);
    }
}

inline void chkpoleAtomI(int i, const MDQCArray<real> x, const MDQCArray<real> y, const MDQCArray<real> z, const MDQCArray<LocalFrame> polaxe,
                         const MDQCArray<int> xaxis, MDQCArray<int> yaxis, const MDQCArray<int> zaxis, MDQCArray2D<real,maxpole> pole)
{
    int k;
    int ia,ib,ic,id;
    real xad,yad,zad;
    real xbd,ybd,zbd;
    real xcd,ycd,zcd;
    real c1,c2,c3,vol;
    bool check;

    check = true;
    if (polaxe[i] != LocalFrame::ZthenX) check = false;
    if (yaxis[i] == 0) check = false;
    if (check) {
        k = yaxis[i];
        ia = i;
        ib = zaxis[i]-1;
        ic = xaxis[i]-1;
        id = std::abs(k)-1;

        // compute the signed parallelpiped volume at chiral site
        xad = x[ia] - x[id];
        yad = y[ia] - y[id];
        zad = z[ia] - z[id];
        xbd = x[ib] - x[id];
        ybd = y[ib] - y[id];
        zbd = z[ib] - z[id];
        xcd = x[ic] - x[id];
        ycd = y[ic] - y[id];
        zcd = z[ic] - z[id];
        c1 = ybd*zcd - zbd*ycd;
        c2 = ycd*zad - zcd*yad;
        c3 = yad*zbd - zad*ybd;
        vol = xad*c1 + xbd*c2 + xcd*c3;
        // invert the multipole components involving the y-axis
        if ((k<0 and vol>0.) or (k>0 and vol<0.)) {
            yaxis[i] = -k;
            pole[i][2] = -pole[i][2];
            pole[i][5] = -pole[i][5];
            pole[i][7] = -pole[i][7];
            pole[i][9] = -pole[i][9];
            pole[i][11] = -pole[i][11];
        }
    }
}
}
