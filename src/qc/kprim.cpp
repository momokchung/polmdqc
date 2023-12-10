// Author: Moses KJ Chung
// Year:   2023

#include "init.h"
#include "kbasis.h"
#include "kprim.h"
#include "mathUtils.h"
#include "unitsqm.h"
#include <iostream>
#include <vector>

namespace prim
{
//////////////////////////////////////
//                                  //
//  kprim  --  primitive gaussians  //
//                                  //
//////////////////////////////////////

// primN          number of primitives
// primLx         x angular momenum
// primLy         y angular momentum
// primLz         z angular momentum
// primExp        exponential of primitive
// primX          x coordinate
// primY          y coordinate
// primZ          z coordinate
// primScale      scale factor
// primNorm       primitive normalization
// primToBasis    map from primitive index to basis index
// primToShell    map from primitive index to shell index

int primN;
// std::vector<int> primLx;
// std::vector<int> primLy;
// std::vector<int> primLz;
// std::vector<real> primExp;
// std::vector<real> primX;
// std::vector<real> primY;
// std::vector<real> primZ;
// std::vector<real> primScale;
std::vector<real> primNorm;
std::vector<int> primToBasis;
std::vector<int> primToShell;

// primShellN        number of primitive shells
// primShellL        total angular momentum
// primShellIndex    primitive index for primitive shell
// primShellX        x coordinate
// primShellY        y coordinate
// primShellZ        z coordinate
// primShellExp      exponential factor for primitive
// primShellCoeff    coefficient for primitive
// primShellScale    primitive scaling

int primShellN;
std::vector<int> primShellL;
std::vector<int> primShellIndex;
std::vector<real> primShellX;
std::vector<real> primShellY;
std::vector<real> primShellZ;
std::vector<real> primShellExp;
std::vector<real> primShellCoeff;
std::vector<real> primShellScale;

void normalizePrimitive();
void buildPrimMaps();

////////////////////////////////////
//                                //
//  kprim  --  set up primitives  //
//                                //
////////////////////////////////////

void kprim()
{
    // initialize
    primN = 0;
    primShellN = 0;
    primShellL.resize(0);
    primShellIndex.resize(0);
    primShellX.resize(0);
    primShellY.resize(0);
    primShellZ.resize(0);
    primShellExp.resize(0);
    primShellCoeff.resize(0);
    primShellScale.resize(0);

    // loop over shell
    for (int i = 0; i < basis::cShellN; ++i)
    {
        int l = basis::cShellL[i];
        real coordx = basis::cShellX[i];
        real coordy = basis::cShellY[i];
        real coordz = basis::cShellZ[i];
        int contractionN = basis::cShellContraction[i];
        auto& coeff = basis::cShellContractionCoeff[i];
        auto& exp = basis::cShellPrimExp[i];
        real scale = basis::cShellScale[i];
        for (int j = 0; j < contractionN; ++j)
        {
            primShellIndex.push_back(primN);
            primN += basis::lToN(l);
            primShellN += 1;
            primShellL.push_back(l);
            primShellX.push_back(coordx);
            primShellY.push_back(coordy);
            primShellZ.push_back(coordz);
            primShellExp.push_back(exp[j]);
            primShellCoeff.push_back(coeff[j]);
            primShellScale.push_back(scale);
        }
    }
    // // print to debug
    // for (int i = 0; i < primShellN; ++i)
    // {
    //     printf("%4i", primShellIndex[i]);
    // }
    // std::cout << std::endl;

    normalizePrimitive();
    // // print to debug
    // for (int i = 0; i < primNorm.size(); ++i)
    // {
    //     printf("%10.5f", primNorm[i]);
    // }
    // std::cout << std::endl;

    buildPrimMaps();
    // // print to debug
    // for (int i = 0; i < primToBasis.size(); ++i)
    // {
    //     printf("%3i", primToBasis[i]);
    // }
    // std::cout << std::endl;
    // for (int i = 0; i < primToShell.size(); ++i)
    // {
    //     printf("%3i", primToShell[i]);
    // }
    // std::cout << std::endl;
}

///////////////////////////////////////////////////
//                                               //
//  normalizePrimitive  --  normalize primitive  //
//                                               //
///////////////////////////////////////////////////

// primitive cartesian gaussian normalization Eq. 2.11 (Fermann, J. T. & Valeev, E. F. Fundamentals of Molecular Integrals Evaluation. (2020).)
void normalizePrimitive()
{
    // primNorm has length primN
    primNorm.resize(0);
    const auto& DF = mathUtils::doubleFactorial;
    auto& partitionL = basis::partitionAngularMomentum;
    real pre1 = pow(2. / unitsqm::pi, 0.75);

    // loop over primitive shell
    for (int i = 0; i < primShellN; ++i)
    {
        int l = primShellL[i];
        real exp = primShellExp[i];
        real pre2 = pre1 * pow(2., l) * pow(exp, l/2. + 0.75);
        for (auto& part: partitionL[l])
        {
            int lx = part[0];
            int ly = part[1];
            int lz = part[2];
            real norm = pre2 / sqrt(DF(2 * lx - 1) * DF(2 * ly - 1) * DF(2 * lz - 1));
            primNorm.push_back(norm);
        }
    }
}

/////////////////////////////////////////////////////////////////////
//                                                                 //
//  buildPrimMaps  --  build prim to basis and prim to shell maps  //
//                                                                 //
/////////////////////////////////////////////////////////////////////

void buildPrimMaps()
{
    // primToBasis has length primN
    primToBasis.resize(0);

    // primToShell has length primShellN
    primToShell.resize(0);
    int basisN = 0;
    auto& partitionL = basis::partitionAngularMomentum;

    // loop over shell
    for (int i = 0; i < basis::cShellN; ++i)
    {
        int l = basis::cShellL[i];
        int contractionN = basis::cShellContraction[i];
        // loop over primitive shell
        for (int j = 0; j < contractionN; ++j)
        {
            int counter = 0;
            // loop over angular momenutum
            for (auto& part: partitionL[l])
            {
                primToBasis.push_back(basisN + counter);
                counter += 1;
            }
            primToShell.push_back(i);
        }
        basisN += basis::lToN(l);
    }
}
}
