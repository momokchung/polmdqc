// Author: Moses KJ Chung
// Year:   2024

#include "alphamol.h"
#include "alphamol1.h"
#include "alphmol.h"
#include "dlauny2.h"

namespace polmdqc
{
////////////////////////////////////////////////////////////
//                                                        //
//  alphamol1  --  single thread surface area and volume  //
//                                                        //
////////////////////////////////////////////////////////////

// "alphamol1" computes volume, surface area, mean, and
// gaussian curvature using a single thread
//
// literature reference:
//
// P. Koehl, A. Akopyan, and H. Edelsbrunner, "Computing the Volume,
// Surface Area, Mean, and Gaussian Curvatures of Molecules and Their
// Derivatives", Journal of Chemical Information and Modeling,
// 63, 973-985, (2023).
//
// github reference:
//
// https://github.com/pkoehl/AlphaMol

void alphamol1(bool deriv)
{
    int natoms = alfatoms.size();
    int nfudge = 8;
    real* surfthd = new real[natoms+nfudge];
    real* volthd = new real[natoms+nfudge];
    real* meanthd = new real[natoms+nfudge];
    real* gaussthd = new real[natoms+nfudge];
    real* dsurfthd = new real[3*(natoms+nfudge)];
    real* dvolthd = new real[3*(natoms+nfudge)];
    real* dmeanthd = new real[3*(natoms+nfudge)];
    real* dgaussthd = new real[3*(natoms+nfudge)];

    alphamol(natoms, &(alfatoms[0]), surfthd, volthd, meanthd, gaussthd,
        dsurfthd, dvolthd, dmeanthd, dgaussthd, deriv);

    wsurf = 0;
    wvol = 0;
    wmean = 0;
    wgauss = 0;
    for(int i = 0; i < natoms; i++) {
        wsurf += surfthd[i];
        wvol += volthd[i];
        wmean += meanthd[i];
        wgauss += gaussthd[i];
    }

    for (int i = 0; i < natoms; i++) {
        surf[alfatoms[i].index] = surfthd[i];
        vol[alfatoms[i].index] = volthd[i];
        mean[alfatoms[i].index] = meanthd[i];
        gauss[alfatoms[i].index] = gaussthd[i];
    }

    if (deriv) {
        for (int i = 0; i < natoms; i++) {
            for (int j = 0; j < 3; j++) {
                dsurf[3*alfatoms[i].index + j] = dsurfthd[3*i+j];
                dvol[3*alfatoms[i].index + j] = dvolthd[3*i+j];
                dmean[3*alfatoms[i].index + j] = dmeanthd[3*i+j];
                dgauss[3*alfatoms[i].index + j] = dgaussthd[3*i+j];
            }
        }
    }

    delete[] surfthd;
    delete[] volthd;
    delete[] meanthd;
    delete[] gaussthd;
    delete[] dsurfthd;
    delete[] dvolthd;
    delete[] dmeanthd;
    delete[] dgaussthd;
}
}
