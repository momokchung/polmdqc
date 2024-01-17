// Author: Moses KJ Chung
// Year:   2023

#include "atoms.h"
#include "bound.h"
#include "boxes.h"
#include "chunks.h"
#include "ewald.h"
#include "fatal.h"
#include "fft.h"
#include "gettext.h"
#include "getword.h"
#include "inform.h"
#include "keys.h"
#include "kewald.h"
#include "lattice.h"
#include "mdqclimits.h"
#include "openmp.h"
#include "pme.h"
#include "potent.h"
#include "upcase.h"
#include <algorithm>
#include <cmath>
#include <sstream>

namespace polmdqc
{
/////////////////////////////////////////////////////
//                                                 //
//  kewald  --  setup for particle mesh Ewald sum  //
//                                                 //
/////////////////////////////////////////////////////

// "kewald" assigns particle mesh Ewald parameters and options
// for a periodic system

void ewaldcof(real& alpha, real& cutoff);
void extent(real& rmax);

void kewald()
{
    constexpr int maxpower = 63;
    constexpr int maxfft = 864;
    int i,k,next;
    int nbig,minfft;
    int iefft1,idfft1;
    int iefft2,idfft2;
    int iefft3,idfft3;
    real delta,rmax;
    real edens,ddens;
    real size,slope;
    real fft1,fft2,fft3;
    std::string keyword;
    std::string record;
    std::string string;
    std::istringstream iss;

    // PME grid size must be even with factors of only 2, 3 and 5
    int multi[maxpower] = {
          2,   4,   6,   8,  10,  12,  16,  18,  20,
         24,  30,  32,  36,  40,  48,  50,  54,  60,
         64,  72,  80,  90,  96, 100, 108, 120, 128,
        144, 150, 160, 162, 180, 192, 200, 216, 240,
        250, 256, 270, 288, 300, 320, 324, 360, 384,
        400, 432, 450, 480, 486, 500, 512, 540, 576,
        600, 640, 648, 720, 750, 768, 800, 810, 864
    };

    // return if Ewald summation is not being used
    if (!use_ewald and !use_dewald) return;

    // set default values for Ewald options and parameters
    ffttyp = "FFTPACK";
    if (nthread > 1) ffttyp = "FFTW";
    boundary = "TINFOIL";
    bseorder = 5;
    bsporder = 5;
    bsdorder = 4;
    edens = 1.2;
    ddens = 0.8;
    aeewald = 0.4;
    apewald = 0.4;
    adewald = 0.4;
    minfft = 16;

    // estimate optimal values for the Ewald coefficient
    if (use_ewald) ewaldcof(aeewald,ewaldcut);
    if (use_dewald) ewaldcof(adewald,dewaldcut);
    if (use_ewald and use_polar) apewald = aeewald;

    // modify Ewald coefficient for small unitcell dimensions
    if (use_polar and use_bounds) {
        size = std::min({xbox,ybox,zbox});
        if (size < 6.) {
            slope = (1.-apewald) / 2.;
            apewald = apewald + slope*(6.-size);
            minfft = 64;
            if (verbose) {
                printf("\n KEWALD  --  Warning, PME Grid Expanded due to Small Cell Size\n");
            }
        }
    }

    // set the system extent for nonperiodic Ewald summation
    if (!use_bounds) {
        extent(rmax);
        xbox = 2. * (rmax+std::max(ewaldcut,dewaldcut));
        ybox = xbox;
        zbox = xbox;
        alphaA = 90.;
        betaA = 90.;
        gammaA = 90.;
        orthogonal = true;
        lattice();
        boundary = "NONE";
        edens = 0.7;
        ddens = 0.7;
    }

    // set defaults for electrostatic and dispersion grid sizes
    nefft1 = 0;
    nefft2 = 0;
    nefft3 = 0;
    ndfft1 = 0;
    ndfft2 = 0;
    ndfft3 = 0;

    // get default grid counts from periodic system dimensions
    delta = 1.0e-8;
    iefft1 = static_cast<int>(xbox*edens-delta) + 1;
    iefft2 = static_cast<int>(ybox*edens-delta) + 1;
    iefft3 = static_cast<int>(zbox*edens-delta) + 1;
    idfft1 = static_cast<int>(xbox*ddens-delta) + 1;
    idfft2 = static_cast<int>(ybox*ddens-delta) + 1;
    idfft3 = static_cast<int>(zbox*ddens-delta) + 1;

    // search keywords for Ewald summation commands
    for (int i = 0; i < nkey; i++) {
        record = keyline[i];
        next = 0;
        upcase(record);
        gettext(record,keyword,next);
        string = record.substr(next);
        iss.clear();
        iss.str(string);
        if (keyword == "FFT-PACKAGE") {
            getword(record,ffttyp,next);
        }
        else if (keyword == "EWALD-ALPHA") {
            iss >> aeewald;
        }
        else if (keyword == "PEWALD-ALPHA") {
            iss >> apewald;
        }
        else if (keyword == "DEWALD-ALPHA") {
            iss >> adewald;
        }
        else if (keyword == "EWALD-BOUNDARY") {
            boundary = "VACUUM";
        }
        else if (keyword == "PME-GRID") {
            fft1 = 0.;
            fft2 = 0.;
            fft3 = 0.;
            iss >> fft1 >> fft2 >> fft3;
            iefft1 = static_cast<int>(std::round(fft1));
            iefft2 = static_cast<int>(std::round(fft2));
            iefft3 = static_cast<int>(std::round(fft3));
            if (iefft2 == 0) iefft2 = iefft1;
            if (iefft3 == 0) iefft3 = iefft1;
        }
        else if (keyword == "DPME-GRID") {
            fft1 = 0.;
            fft2 = 0.;
            fft3 = 0.;
            iss >> fft1 >> fft2 >> fft3;;
            idfft1 = static_cast<int>(std::round(fft1));;
            idfft2 = static_cast<int>(std::round(fft2));;
            idfft3 = static_cast<int>(std::round(fft3));;
            if (idfft2 == 0) idfft2 = idfft1;
            if (idfft3 == 0) idfft3 = idfft1;
        }
        else if (keyword == "PME-ORDER") {
            iss >> bseorder;
        }
        else if (keyword == "PPME-ORDER") {
            iss >> bsporder;
        }
        else if (keyword == "DPME-ORDER") {
            iss >> bsdorder;
        }
    }

    // determine electrostatic grid size from allowed values
    if (use_ewald) {
        nefft1 = maxfft;
        nefft2 = maxfft;
        nefft3 = maxfft;
        for (int i = maxpower-1; i >= 0; i--) {
            k = multi[i];
            if (k <= maxfft) {
                if (k >= iefft1) nefft1 = k;
                if (k >= iefft2) nefft2 = k;
                if (k >= iefft3) nefft3 = k;
            }
        }
        if (nefft1 < minfft) nefft1 = minfft;
        if (nefft2 < minfft) nefft2 = minfft;
        if (nefft3 < minfft) nefft3 = minfft;
    }

    // determine dispersion grid size from allowed values
    if (use_dewald) {
        ndfft1 = maxfft;
        ndfft2 = maxfft;
        ndfft3 = maxfft;
        for (int i = maxpower-1; i >= 0; i--) {
            k = multi[i];
            if (k <= maxfft) {
                if (k >= idfft1) ndfft1 = k;
                if (k >= idfft2) ndfft2 = k;
                if (k >= idfft3) ndfft3 = k;
            }
        }
        if (ndfft1 < minfft) ndfft1 = minfft;
        if (ndfft2 < minfft) ndfft2 = minfft;
        if (ndfft3 < minfft) ndfft3 = minfft;
    }

    // check the particle mesh Ewald grid dimensions
    nbig = std::max({nefft1,nefft2,nefft3,ndfft1,ndfft2,ndfft3});
    if (nbig > maxfft) {
        printf("\n KEWALD  --  PME Grid Size Too Large; Increase MAXFFT\n");
        fatal();
    }
    if (use_ewald and (nefft1<iefft1 or nefft2<iefft2 or nefft3<iefft3)) {
        printf("\n KEWALD  --  Warning, Small ElectrostaticPME Grid Size\n");
    }
    if (use_dewald and (ndfft1<idfft1 or ndfft2<idfft2 or ndfft3<idfft3)) {
        printf("\n KEWALD  --  Warning, Small DispersionPME Grid Size\n");
    }

    // set maximum sizes for PME grid and B-spline order
    nfft1 = std::max(nefft1,ndfft1);
    nfft2 = std::max(nefft2,ndfft2);
    nfft3 = std::max(nefft3,ndfft3);
    bsorder = std::max({bseorder,bsporder,bsdorder});

    // perform dynamic allocation of some global arrays
    if (bsmod1.size() != 0) bsmod1.resize(0);
    if (bsmod2.size() != 0) bsmod2.resize(0);
    if (bsmod3.size() != 0) bsmod3.resize(0);
    if (bsbuild.size() != 0) bsbuild.resize(0);
    if (thetai1.size() != 0) thetai1.resize(0);
    if (thetai2.size() != 0) thetai2.resize(0);
    if (thetai3.size() != 0) thetai3.resize(0);
    if (pmetable.size() != 0) pmetable.resize(0);
    bsmod1.resize(nfft1);
    bsmod2.resize(nfft2);
    bsmod3.resize(nfft3);
    bsbuild.resize(bsorder, std::vector<real>(bsorder));
    thetai1.resize(n, std::vector<std::vector<real>>(bsorder, std::vector<real>(4)));
    thetai2.resize(n, std::vector<std::vector<real>>(bsorder, std::vector<real>(4)));
    thetai3.resize(n, std::vector<std::vector<real>>(bsorder, std::vector<real>(4)));
    pmetable.resize(6*nthread, std::vector<int>(n));

    // print a message listing some of the Ewald parameters
    if (verbose) {
        printf("\n Particle Mesh Ewald Parameters :\n\n");
        printf("     Type                Ewald Alpha    Grid Dimensions    Spline Order\n\n");
        if (use_ewald) {
            printf("   Electrostatics         %8.4f     %5d%5d%5d       %5d\n", aeewald,nefft1,nefft2,nefft3,bseorder);
            if (use_polar) {
                printf("   Polarization           %8.4f     %5d%5d%5d       %5d\n", apewald,nefft1,nefft2,nefft3,bsporder);
            }
        }
        if (use_dewald) {
            printf("   Dispersion             %8.4f     %5d%5d%5d       %5d\n", adewald,ndfft1,ndfft2,ndfft3,bsdorder);
        }
    }
}

/////////////////////////////////////////////////////
//                                                 //
//  ewaldcof  --  estimation of Ewald coefficient  //
//                                                 //
/////////////////////////////////////////////////////

// "ewaldcof" finds an Ewald coefficient such that all terms
// beyond the specified cutoff distance will have a value less
// than a specified tolerance


void ewaldcof(real& alpha, real& cutoff)
{
    int i,k;
    real x,xlo,xhi,y,ratio,eps;

    // set tolerance value; use of 1.0d-8 over 1.0d-6 gives
    // larger Ewald coefficients to ensure gradient continuity
    eps = 1.0e-8;

    // get approximate value from cutoff and tolerance
    ratio = eps + 1.;
    x = 0.5;
    i = 0;
    while (ratio >= eps) {
        i++;
        x = 2. * x;
        y = x * cutoff;
        ratio = std::erfc(y) / cutoff;
    }

    // use a binary search to refine the coefficient
    k = i + 60;
    xlo = 0.;
    xhi = x;
    for (int i = 0; i < k; i++) {
        x = (xlo+xhi) / 2.;
        y = x * cutoff;
        ratio = std::erfc(y) / cutoff;
        if (ratio >= eps) xlo = x;
        else xhi = x;
    }
    alpha = x;
}

/////////////////////////////////////////////////////
//                                                 //
//  extent  --  find maximum interatomic distance  //
//                                                 //
/////////////////////////////////////////////////////

// "extent" finds the largest interatomic distance in a system


void extent(real& rmax)
{
    int i,k;
    real xi,yi,zi,xk,yk,zk,r2;

    // search all atom pairs to find the largest distance
    rmax = 0.;
    for (int i = 0; i < n-1; i++) {
        xi = x[i];
        yi = y[i];
        zi = z[i];
        for (int k = i+1; k < n; k++) {
            xk = x[k];
            yk = y[k];
            zk = z[k];
            r2 = std::pow((xk-xi),2) + std::pow((yk-yi),2) + std::pow((zk-zi),2);
            rmax = std::max(r2,rmax);
        }
    }
    rmax = std::sqrt(rmax);
}
}
