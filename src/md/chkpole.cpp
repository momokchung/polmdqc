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

void chkpole()
{
    int i,k;
    int ia,ib,ic,id;
    double xad,yad,zad;
    double xbd,ybd,zbd;
    double xcd,ycd,zcd;
    double c1,c2,c3,vol;
    bool dopol,dorep;
    bool check;

    // loop over multipoles and test for chirality inversion
    for (int i = 0; i < n; i++) {
        dopol = false;
        dorep = false;
        if (pollist.size() != 0) {
            if (pollist[0] != -1) dopol = true;
        }
        if (replist.size() != 0) {
            if (replist[i] != -1) dorep = true;
        }
        if (dopol or dorep) {
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
                    if (dopol) {
                        pole[i][2] = -pole[i][2];
                        pole[i][5] = -pole[i][5];
                        pole[i][7] = -pole[i][7];
                        pole[i][9] = -pole[i][9];
                        pole[i][11] = -pole[i][11];
                    }
                    if (dorep) {
                        repole[i][2] = -repole[i][2];
                        repole[i][5] = -repole[i][5];
                        repole[i][6] = -repole[i][6];
                        repole[i][9] = -repole[i][9];
                        repole[i][11] = -repole[i][11];
                    }
                }
            }
        }
    }
}
}
