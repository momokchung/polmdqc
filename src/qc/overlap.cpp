// Author: Moses KJ Chung
// Year:   2023

#include "kbasis.h"
#include "kgbs.h"
#include "kprim.h"
#include "mathUtils.h"
#include "overlap.h"
#include "print.h"
#include "unitsqm.h"
#include <iostream>
#include <vector>

namespace overlap
{
///////////////////////////////////
//                               //
//  overlap  --  overlap matrix  //
//                               //
///////////////////////////////////

// cartS    cartesian overlap matrix
// sphS     spherical overlap matrix

std::vector<std::vector<real>> cartS;
std::vector<std::vector<real>> sphS;

////////////////////////////////////////////////////////////////////////
//                                                                    //
//  overlapOS  --  construct overlap matrix using Obara Saika method  //
//                                                                    //
////////////////////////////////////////////////////////////////////////

void overlapOS()
{
    int basisN = basis::basisN;
    int maxL = basis::cShellLMax;
    int shellN = prim::primShellN;
    auto& partitionL = basis::partitionAngularMomentum;

    // initialize xS, yS, and zS matrix
    std::vector<std::vector<real>> xS(maxL + 1, std::vector <real>(maxL + 1));
    std::vector<std::vector<real>> yS(maxL + 1, std::vector <real>(maxL + 1));
    std::vector<std::vector<real>> zS(maxL + 1, std::vector <real>(maxL + 1));

    // allocate and initialize overlap matrix
    cartS.resize(0);
    cartS.resize(basisN, std::vector<real>(basisN, 0.));

    // outer loop over primitive shell
    for (int i = 0; i < shellN; ++i)
    {
        int iL = prim::primShellL[i];
        int iIndex = prim::primShellIndex[i];
        int iShell = prim::primToShell[i];
        real ai = prim::primShellExp[i];
        real coeffi = prim::primShellCoeff[i];
        real coordxi = prim::primShellX[i];
        real coordyi = prim::primShellY[i];
        real coordzi = prim::primShellZ[i];
        
        // inner loop over primitive shell
        for (int j = 0; j <= i; ++j)
        {
            int jL = prim::primShellL[j];
            int jIndex = prim::primShellIndex[j];
            int jShell = prim::primToShell[j];
            real aj = prim::primShellExp[j];
            real coeffj = prim::primShellCoeff[j];
            real coordxj = prim::primShellX[j];
            real coordyj = prim::primShellY[j];
            real coordzj = prim::primShellZ[j];

            real coordx = coordxj - coordxi;
            real coordy = coordyj - coordyi;
            real coordz = coordzj - coordzi;
            real r2 = coordx * coordx + coordy * coordy + coordz * coordz;

            real aij = ai * aj;
            real aP = ai + aj;
            real aP2 = 2. * aP;
            real xP = (ai*coordxi + aj*coordxj)/aP;
            real yP = (ai*coordyi + aj*coordyj)/aP;
            real zP = (ai*coordzi + aj*coordzj)/aP;
            real xPI = xP - coordxi;
            real yPI = yP - coordyi;
            real zPI = zP - coordzi;
            real xPJ = xP - coordxj;
            real yPJ = yP - coordyj;
            real zPJ = zP - coordzj;
            real kp = exp(-aij/aP*r2);
            real pre = kp * coeffi * coeffj * pow(unitsqm::pi/aP,1.5);

            // initialize xS, yS, zS matrix
            for (int k = 0; k <= iL; ++k)
            {
                for (int l = 0; l <= jL; ++l)
                {
                    xS[k][l] = 0.;
                    yS[k][l] = 0.;
                    zS[k][l] = 0.;
                }
            }
            xS[0][0] = 1.;
            yS[0][0] = 1.;
            zS[0][0] = 1.;

            // build xS, yS, zS matrix
            for (int k = 0; k <= iL; ++k)
            {
                recursionOA(k, 0, xPI, aP2, xS);
                recursionOA(k, 0, yPI, aP2, yS);
                recursionOA(k, 0, zPI, aP2, zS);
                for (int l = 0; l <= jL; ++l)
                {
                    recursionOB(k, l, xPJ, aP2, xS);
                    recursionOB(k, l, yPJ, aP2, yS);
                    recursionOB(k, l, zPJ, aP2, zS);
                }
            }

            // increment primitive overlap to basis overlap
            int iPrim = iIndex;
            for (auto& parti: partitionL[iL])
            {
                int lxi = parti[0];
                int lyi = parti[1];
                int lzi = parti[2];
                int basisNi = prim::primToBasis[iPrim];
                int jPrim = jIndex;
                for (auto& partj: partitionL[jL])
                {
                    int lxj = partj[0];
                    int lyj = partj[1];
                    int lzj = partj[2];
                    int basisNj = prim::primToBasis[jPrim];
                    real sIntermediate = xS[lxi][lxj] * yS[lyi][lyj] * zS[lzi][lzj];
                    if ((i != j) and (iShell == jShell))
                        sIntermediate *= 2.;
                    cartS[basisNi][basisNj] += pre * prim::primNorm[iPrim] * prim::primNorm[jPrim] * sIntermediate;
                    jPrim += 1;
                }
                iPrim += 1;
            }
        }
    }

    // loop over basis for normalization
    for (int i = 0; i < basisN; ++i)
    {
        real normi = basis::basisNorm[i];
        for (int j = 0; j <= i; ++j)
        {
            cartS[i][j] *= normi * basis::basisNorm[j];
        }
    }

    // symmetrize cartesian overlap
    mathUtils::symmetrize(cartS);

    // construct spherical overlap
    if (gbs::basisType == gbs::BasisType::spherical)
    {
        int sphBasisN = basis::sphBasisN;
        sphS.resize(0);
        sphS.resize(sphBasisN, std::vector<real>(sphBasisN, 0.));

        // loop over cartesian basis, add to appropriate spherical basis
        for (int i = 0; i < basisN; ++i)
        {
            auto& iContraction = basis::cartSphContraction[i];
            auto& iCoeff = basis::cartSphCoeff[i];
            int iN = iContraction.size();

            for (int j = 0; j < basisN; ++j)
            {
                auto& jContraction = basis::cartSphContraction[j];
                auto& jCoeff = basis::cartSphCoeff[j];
                int jN = jContraction.size();

                real s = cartS[i][j];

                for (int ii = 0; ii < iN; ++ii)
                {
                    int cii = iContraction[ii];
                    real coefii = iCoeff[ii];

                    for (int jj = 0; jj < jN; ++jj)
                    {
                        int cjj = jContraction[jj];
                        real coefjj = jCoeff[jj];

                        sphS[cii][cjj] += coefii * coefjj * s;
                    }
                }
            }
        }
    }

    // // print to debug
    // if (gbs::basisType == gbs::BasisType::cartesian)
    // {
    //     print::printMatrix(cartS, "Cartesian Overlap Integral");
    // }
    // else if (gbs::basisType == gbs::BasisType::spherical)
    // {
    //     print::printMatrix(sphS, "Spherical Overlap Integral");
    // }
}
}
