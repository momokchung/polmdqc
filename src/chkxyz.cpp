////////////////////////////////////////////////////////
//                                                    //
//  chkxyz.cpp  --  check for coincident coordinates  //
//                                                    //
////////////////////////////////////////////////////////

// "chkxyz" finds any pairs of atoms with identical Cartesian
// coordinates, and prints a warning message


#include "atoms.h"
#include "chkxyz.h"

void chkxyz(bool& clash)
{
    // initialize distance tolerance and atom collision flag
    double eps = 0.000001;
    clash = false;
    bool header = true;

    // loop over atom pairs testing for identical coordinates
    for (int i = 0; i < n-1; i++) {
        double xi = x[i];
        double yi = y[i];
        double zi = z[i];
        for (int k = i+1; k < n; k++) {
            double xr = x[k] - xi;
            double yr = y[k] - yi;
            double zr = z[k] - zi;
            double r2 = xr*xr + yr*yr + zr*zr;
            if (r2 < eps) {
                clash = true;
                if (header) {
                    header = false;
                    printf("\n");
                }
                printf("\n CHKXYZ  --  Warning, Atoms%6d and%6d have Identical Coordinates\n", i+1, k+1);
            }
        }
    }
}
