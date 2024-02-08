// Author: Moses KJ Chung
// Year:   2023

#include "boys.h"
#include "katomsqm.h"
#include "kbasis.h"
#include "kgbs.h"
#include "kprim.h"
#include "mathUtils.h"
#include "nuclear.h"
#include "print.h"
#include "unitsqm.h"
#include <iostream>
#include <vector>

namespace nuclear
{
///////////////////////////////////////////////////////
//                                                   //
//  nuclear  --  nuclear electron attraction matrix  //
//                                                   //
///////////////////////////////////////////////////////

// cartNE    cartesian nuclear electron attraction matrix
// sphNE     spherical nuclear electron attraction matrix

std::vector<std::vector<real>> cartNE;
std::vector<std::vector<real>> sphNE;

////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                        //
//  nuclearOS  --  construct nuclear electron attraction matrix using Obara Saika method  //
//                                                                                        //
////////////////////////////////////////////////////////////////////////////////////////////

void nuclearOS()
{
    // get atom information
    int atomN = atoms::n;

    // get basis and primitive information
    int basisN = basis::basisN;
    int maxL = basis::cShellLMax;
    int shellN = prim::primShellN;
    auto& partitionL = basis::partitionAngularMomentum;

    // initialize intermediate vectors
    std::vector<real> poly(2 * maxL + 1);
    std::vector<real> polyxy(2 * maxL + 1);
    std::vector<real> boysPoly(2 * maxL + 1);

    // initialize xNe, yNe, and zNe matrix
    int maxTOrder = 2 * maxL + 1;
    std::vector<std::vector<std::vector<real>>> xNe(maxL + 1, std::vector<std::vector<real>>(maxL + 1, std::vector<real>(maxTOrder)));
    std::vector<std::vector<std::vector<real>>> yNe(maxL + 1, std::vector<std::vector<real>>(maxL + 1, std::vector<real>(maxTOrder)));
    std::vector<std::vector<std::vector<real>>> zNe(maxL + 1, std::vector<std::vector<real>>(maxL + 1, std::vector<real>(maxTOrder)));

    // allocate and initialize nuclear electron attraction matrix
    cartNE.resize(0);
    cartNE.resize(basisN, std::vector<real>(basisN, 0.));

    // outer loop over primitive shell
    for (int i = 0; i < shellN; i++)
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
        for (int j = 0; j <= i; j++)
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
            real pre = kp * coeffi * coeffj * 2. * unitsqm::pi / aP;

            // loop over nuclei
            for (int k = 0; k < atomN; k++)
            {
                int charge = atoms::core[k];
                real prefactor = -charge * pre;
                real coordxk = atoms::coordx[k];
                real coordyk = atoms::coordy[k];
                real coordzk = atoms::coordz[k];
                
                real xPK = xP - coordxk;
                real yPK = yP - coordyk;
                real zPK = zP - coordzk;

                // initialize xNe, yNe, zNe matrix
                for (int l = 0; l <= iL; l++)
                {
                    for (int m = 0; m <= jL; m++)
                    {
                        int torder = l + m;
                        for (int n = 0; n <= torder; n++)
                        {
                            xNe[l][m][n] = 0.;
                            yNe[l][m][n] = 0.;
                            zNe[l][m][n] = 0.;
                        }
                    }
                }
                xNe[0][0][0] = 1.;
                yNe[0][0][0] = 1.;
                zNe[0][0][0] = 1.;

                // build xNe, yNe, zNe matrix
                for (int l = 0; l <= iL; l++)
                {
                    recursionNeA(l, 0, xPI, xPK, aP2, xNe);
                    recursionNeA(l, 0, yPI, yPK, aP2, yNe);
                    recursionNeA(l, 0, zPI, zPK, aP2, zNe);
                    for (int m = 0; m <= jL; m++)
                    {
                        recursionNeB(l, m, xPJ, xPK, aP2, xNe);
                        recursionNeB(l, m, yPJ, yPK, aP2, yNe);
                        recursionNeB(l, m, zPJ, zPK, aP2, zNe);
                    }
                }

                // get boys integral
                real T = aP * (xPK * xPK + yPK * yPK + zPK * zPK);
                boys::boysIntegralPoly(T, iL + jL, boysPoly);
                
                // increment primitive nuclear attraction to basis nuclear attraction
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

                        int lx = lxi + lxj;
                        int ly = lyi + lyj;
                        int lz = lzi + lzj;

                        std::vector<real>& polyx = xNe[lxi][lxj];
                        std::vector<real>& polyy = yNe[lyi][lyj];
                        std::vector<real>& polyz = zNe[lzi][lzj];

                        // zero out polynomials
                        mathUtils::zero(lx + ly + 1, polyxy);
                        mathUtils::zero(lx + ly + lz + 1, poly);

                        // multiply polynomials
                        mathUtils::multpoly(lx + 1, ly + 1, polyx, polyy, polyxy);
                        mathUtils::multpoly(lx + ly + 1, lz + 1, polyxy, polyz, poly);

                        // contract polynomials
                        real neIntermediate = mathUtils::contract(iL + jL + 1, poly, boysPoly);

                        if ((i != j) and (iShell == jShell))
                            neIntermediate *= 2.;
                        cartNE[basisNi][basisNj] += prefactor * prim::primNorm[iPrim] * prim::primNorm[jPrim] * neIntermediate;
                        jPrim += 1;
                    }
                    iPrim += 1;
                }
            }
        }
    }

    // loop over basis for normalization
    for (int i = 0; i < basisN; i++)
    {
        real normi = basis::basisNorm[i];
        for (int j = 0; j <= i; j++)
        {
            cartNE[i][j] *= normi * basis::basisNorm[j];
        }
    }

    // symmetrize nuclear electron attraction
    mathUtils::symmetrize(cartNE);

    // construct spherical nuclear electron attraction
    if (gbs::basisType == gbs::BasisType::spherical)
    {
        int sphBasisN = basis::sphBasisN;
        sphNE.resize(0);
        sphNE.resize(sphBasisN, std::vector<real>(sphBasisN, 0.));

        // loop over cartesian basis, add to appropriate spherical basis
        for (int i = 0; i < basisN; i++)
        {
            auto& iContraction = basis::cartSphContraction[i];
            auto& iCoeff = basis::cartSphCoeff[i];
            int iN = iContraction.size();

            for (int j = 0; j < basisN; j++)
            {
                auto& jContraction = basis::cartSphContraction[j];
                auto& jCoeff = basis::cartSphCoeff[j];
                int jN = jContraction.size();

                real ne = cartNE[i][j];

                for (int ii = 0; ii < iN; ii++)
                {
                    int cii = iContraction[ii];
                    real coefii = iCoeff[ii];

                    for (int jj = 0; jj < jN; jj++)
                    {
                        int cjj = jContraction[jj];
                        real coefjj = jCoeff[jj];

                        sphNE[cii][cjj] += coefii * coefjj * ne;
                    }
                }
            }
        }
    }

    // // print to debug
    // if (gbs::basisType == gbs::BasisType::cartesian)
    // {
    //     print::printMatrix(cartNE, "Cartesian Nuclear Attraction Integral");
    // }
    // else if (gbs::basisType == gbs::BasisType::spherical)
    // {
    //     print::printMatrix(sphNE, "Spherical Nuclear Attraction Integral");
    // }
}
}
