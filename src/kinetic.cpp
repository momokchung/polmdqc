//////////////////////////////////////////////////////
//                                                  //
//  kinetic.cpp  --  compute kinetic energy matrix  //
//                                                  //
//////////////////////////////////////////////////////


#include "kbasis.h"
#include "kgbs.h"
#include "kinetic.h"
#include "kprim.h"
#include "mathUtils.h"
#include "overlap.h"
#include "print.h"
#include "unitsqm.h"
#include <iostream>
#include <vector>

namespace kinetic
{
// cartKE    cartesian kinetic energy matrix
// sphKE     spherical kinetic energy matrix

std::vector<std::vector<real>> cartKE;
std::vector<std::vector<real>> sphKE;


//////////////////////////////////////////////////////////////////////////////////
//                                                                              //
//  void kineticOS  --  compute kinetic energy matrix using Obara Saika method  //
//                                                                              //
//////////////////////////////////////////////////////////////////////////////////

void kineticOS()
{
    int basisN = basis::basisN;
    int maxL = basis::cShellLMax;
    int shellN = prim::primShellN;
    auto& partitionL = basis::partitionAngularMomentum;

    // initialize xS, yS, and zS matrix
    std::vector<std::vector<real>> xS(maxL + 2, std::vector<real>(maxL + 2));
    std::vector<std::vector<real>> yS(maxL + 2, std::vector<real>(maxL + 2));
    std::vector<std::vector<real>> zS(maxL + 2, std::vector<real>(maxL + 2));

    // initialize xKE, yKE, and zKE matrix
    std::vector<std::vector<real>> xKE(maxL + 1, std::vector<real>(maxL + 1));
    std::vector<std::vector<real>> yKE(maxL + 1, std::vector<real>(maxL + 1));
    std::vector<std::vector<real>> zKE(maxL + 1, std::vector<real>(maxL + 1));

    // allocate and initialize kinetic energy matrix
    cartKE.resize(0);
    cartKE.resize(basisN, std::vector<real>(basisN, 0.));

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
            for (int k = 0; k <= (iL + 1); ++k)
            {
                for (int l = 0; l <= (jL +1); ++l)
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
            for (int k = 0; k <= (iL + 1); ++k)
            {
                overlap::recursionOA(k, 0, xPI, aP2, xS);
                overlap::recursionOA(k, 0, yPI, aP2, yS);
                overlap::recursionOA(k, 0, zPI, aP2, zS);
                for (int l = 0; l <= (jL + 1); ++l)
                {
                    overlap::recursionOB(k, l, xPJ, aP2, xS);
                    overlap::recursionOB(k, l, yPJ, aP2, yS);
                    overlap::recursionOB(k, l, zPJ, aP2, zS);
                }
            }
            
            // initialize xKE, yKE, zKE matrix
            for (int k = 0; k <= iL; ++k)
            {
                for (int l = 0; l <= jL; ++l)
                {
                    xKE[k][l] = 0.;
                    yKE[k][l] = 0.;
                    zKE[k][l] = 0.;
                }
            }

            // build xKE, yKE, zKE matrix
            for (int k = 0; k <= iL; ++k)
            {
                
                for (int l = 0; l <= jL; ++l)
                {
                    xKE[k][l] += 2. * ai * aj * xS[k + 1][l + 1];
                    yKE[k][l] += 2. * ai * aj * yS[k + 1][l + 1];
                    zKE[k][l] += 2. * ai * aj * zS[k + 1][l + 1];
                    if (k > 0)
                    {
                        xKE[k][l] -= k * aj * xS[k - 1][l + 1];
                        yKE[k][l] -= k * aj * yS[k - 1][l + 1];
                        zKE[k][l] -= k * aj * zS[k - 1][l + 1];
                    }
                    if (l > 0)
                    {
                        xKE[k][l] -= l * ai * xS[k + 1][l - 1];
                        yKE[k][l] -= l * ai * yS[k + 1][l - 1];
                        zKE[k][l] -= l * ai * zS[k + 1][l - 1];
                    }
                    if (k > 0 and l > 0)
                    {
                        xKE[k][l] += 0.5 * k * l * xS[k - 1][l - 1];
                        yKE[k][l] += 0.5 * k * l * yS[k - 1][l - 1];
                        zKE[k][l] += 0.5 * k * l * zS[k - 1][l - 1];
                    }
                }
            }

            // increment primitive kinetic energy to basis kinetic energy
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
                    real keIntermediatex = xKE[lxi][lxj] * yS[lyi][lyj] * zS[lzi][lzj];
                    real keIntermediatey = xS[lxi][lxj] * yKE[lyi][lyj] * zS[lzi][lzj];
                    real keIntermediatez = xS[lxi][lxj] * yS[lyi][lyj] * zKE[lzi][lzj];
                    real keIntermediate = keIntermediatex + keIntermediatey + keIntermediatez;

                    if ((i != j) && (iShell == jShell))
                        keIntermediate *= 2.;
                    cartKE[basisNi][basisNj] += pre * prim::primNorm[iPrim] * prim::primNorm[jPrim] * keIntermediate;
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
            cartKE[i][j] *= normi * basis::basisNorm[j];
        }
    }

    // symmetrize kinetic energy
    mathUtils::symmetrize(cartKE);

    // construct spherical kinetic energy
    if (gbs::basisType == gbs::BasisType::spherical)
    {
        int sphBasisN = basis::sphBasisN;
        sphKE.resize(0);
        sphKE.resize(sphBasisN, std::vector<real>(sphBasisN, 0.));

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

                real ke = cartKE[i][j];

                for (int ii = 0; ii < iN; ++ii)
                {
                    int cii = iContraction[ii];
                    real coefii = iCoeff[ii];

                    for (int jj = 0; jj < jN; ++jj)
                    {
                        int cjj = jContraction[jj];
                        real coefjj = jCoeff[jj];

                        sphKE[cii][cjj] += coefii * coefjj * ke;
                    }
                }
            }
        }
    }

    // // print to debug
    // if (gbs::basisType == gbs::BasisType::cartesian)
    // {
    //     print::printMatrix(cartKE, "Cartesian Kinetic-Energy Integral");
    // }
    // else if (gbs::basisType == gbs::BasisType::spherical)
    // {
    //     print::printMatrix(sphKE, "Spherical Kinetic-Energy Integral");
    // }
}
}
