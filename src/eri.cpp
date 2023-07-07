//////////////////////////////////////////////
//                                          //
//  eri.cpp  --  electron repulsion matrix  //
//                                          //
//////////////////////////////////////////////


#include "boys.h"
#include "config.h"
#include "kbasis.h"
#include "kgbs.h"
#include "kprim.h"
#include "mathUtils.h"
#include "eri.h"
#include "units.h"
#include <iostream>
#include <vector>

namespace eri
{
// cartERI    cartesian electron repulsion matrix
// sphERI     spherical electron repulsion matrix

std::vector<real> cartERI;
std::vector<real> sphERI;
std::vector<std::vector<real>> pairPrimPx;
std::vector<std::vector<real>> pairPrimPy;
std::vector<std::vector<real>> pairPrimPz;
std::vector<std::vector<real>> pairPrimZ;
std::vector<std::vector<real>> pairPrimK;


//////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                      //
//  void inCoreEriOS  --  construct electron repulsion matrix using Obara Saika method for in core SCF  //
//                                                                                                      //
//////////////////////////////////////////////////////////////////////////////////////////////////////////

void inCoreEriOS()
{
    int basisN = basis::basisN;
    int maxL = basis::cShellLMax;
    int shellN = prim::primShellN;
    auto& partitionL = basis::partitionAngularMomentum;

    // initialize ... matrices

    // check if we have enough RAM
    int memory = config::memory * 0.9;
    int basisN4 = basisN * basisN * basisN * basisN;
    int eriMemory = basisN4 * 8 / 1000000;
    if ((eriMemory) > memory)
    {
        std::cerr << eriMemory << " MB is needed for ERI RAM storage but only "
                  << memory << " MB was provided" << std::endl;
        std::exit(1);
    }

    // allocate and initialize primitive pair matrices
    pairPrimPx.resize(shellN, std::vector<real>(shellN, 0.));
    pairPrimPy.resize(shellN, std::vector<real>(shellN, 0.));
    pairPrimPz.resize(shellN, std::vector<real>(shellN, 0.));
    pairPrimZ.resize(shellN, std::vector<real>(shellN, 0.));
    pairPrimK.resize(shellN, std::vector<real>(shellN, 0.));

    // allocate and initialize ERI matrix
    int pairsN = basisN * (basisN + 1) / 2;
    int uniqueM = pairsN * (pairsN + 1) / 2;
    cartERI.resize(uniqueM, 0.);

    // outer loop over primitive shells
    for (int i = 0; i < shellN; ++i)
    {
        real ai = prim::primShellExp[i];
        real coordxi = prim::primShellX[i];
        real coordyi = prim::primShellY[i];
        real coordzi = prim::primShellZ[i];

        // inner loop over primitive shell
        for (int j = 0; j <= i; ++j)
        {
            real aj = prim::primShellExp[j];
            real coordxj = prim::primShellX[j];
            real coordyj = prim::primShellY[j];
            real coordzj = prim::primShellZ[j];

            real coordx = coordxj - coordxi;
            real coordy = coordyj - coordyi;
            real coordz = coordzj - coordzi;
            real r2 = coordx * coordx + coordy * coordy + coordz * coordz;

            real aij = ai * aj;
            real aP = ai + aj;

            real xP = (ai*coordxi + aj*coordxj)/aP;
            real yP = (ai*coordyi + aj*coordyj)/aP;
            real zP = (ai*coordzi + aj*coordzj)/aP;

            real k = sqrt(2.) * pow(units::pi,1.25) / aP * exp(-aij / aP * r2);

            pairPrimPx[i][j] = xP;
            pairPrimPy[i][j] = yP;
            pairPrimPz[i][j] = zP;
            pairPrimZ[i][j] = aP;
            pairPrimK[i][j] = k;
        }
    }

    int primPairN = shellN * (shellN + 1) / 2;

    int i = 0;
    int j = 0;
    // outer loop over pair of primitive shells
    for (int ipps = 0; ipps < primPairN; ++ipps)
    {
        real xP1 = pairPrimPx[i][j];
        real yP1 = pairPrimPy[i][j];
        real zP1 = pairPrimPz[i][j];
        real aP1 = pairPrimZ[i][j];
        real k1 = pairPrimK[i][j];

        // inner loop over pair of primitive shells
        int k = 0;
        int l = 0;
        for (int jpps = 0; jpps <= ipps; ++jpps)
        {
            real xP2 = pairPrimPx[k][l];
            real yP2 = pairPrimPy[k][l];
            real zP2 = pairPrimPz[k][l];
            real aP2 = pairPrimZ[k][l];
            real k2 = pairPrimK[k][l];

            // increment k and l
            k += 1;
            if ((k % shellN) == 0)
            {
                l += 1;
                k = l;
            }
        }

        // increment i and j
        i += 1;
        if ((i % shellN) == 0)
        {
            j += 1;
            i = j;
        }
    }
}

}
