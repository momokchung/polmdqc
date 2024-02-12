// Author: Moses KJ Chung
// Year:   2023

#include "atoms.h"
#include "chkxyz.h"

namespace polmdqc
{
////////////////////////////////////////////////////
//                                                //
//  chkxyz  --  check for coincident coordinates  //
//                                                //
////////////////////////////////////////////////////

// "chkxyz" finds any pairs of atoms with identical Cartesian
// coordinates, and prints a warning message

void chkxyz(bool& clash)
{
    // initialize distance tolerance and atom collision flag
    real eps = 0.000001;
    clash = false;
    bool header = true;

    // loop over atom pairs testing for identical coordinates
    for (int i = 0; i < n-1; i++) {
        real xi = x[i];
        real yi = y[i];
        real zi = z[i];
        for (int k = i+1; k < n; k++) {
            real xr = x[k] - xi;
            real yr = y[k] - yi;
            real zr = z[k] - zi;
            real r2 = xr*xr + yr*yr + zr*zr;
            if (r2 < eps) {
                clash = true;
                if (header) {
                    header = false;
                    printf("\n");
                }
                printf(" CHKXYZ  --  Warning, Atoms%6d and%6d have Identical Coordinates\n", i+1, k+1);
            }
        }
    }
}
}
