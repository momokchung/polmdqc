/////////////////////////////////
//                             //
//  kbasis.cpp  --  basis set  //
//                             //
/////////////////////////////////


#include "kbasis.h"
#include "katomsqm.h"
#include "kgbs.h"
#include "mathUtils.h"
#include "units.h"
#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

namespace basis
{
// basisN        number of cartesian basis
// basisLx       x angular momentum
// basisLy       y angular momenum
// basisLz       z angular momentum
// basisL        total angular momentum
// basisCoeff    primitive coefficients
// basisExp      primitive exponents
// basisX        x coordinate
// basisY        y coordinate
// basisZ        z coordinate
// basisNorm     cartesian basis normalization

int basisN = 0;
// std::vector<int> basisLx;
// std::vector<int> basisLy;
// std::vector<int> basisLz;
// std::vector<int> basisL;
// std::vector<std::vector<real>> basisCoeff;
// std::vector<std::vector<real>> basisExp;
// std::vector<real> basisX;
// std::vector<real> basisY;
// std::vector<real> basisZ;
std::vector<real> basisNorm;

// sphBasisN             number of spherical basis
// sphContraction        returns list of cartesian basis index for a given spherical basis index
// sphCoeff              returns list of carteisan coefficients for a given spherical basis
// cartSphContraction    returns list of spherical basis index for a given cartesian basis index
// cartSphCoeff          returns list of spherical coefficients for a given cartesian basis

int sphBasisN = 0;
std::vector<std::vector<int>> sphContraction;
std::vector<std::vector<real>> sphCoeff;
std::vector<std::vector<int>> cartSphContraction;
std::vector<std::vector<real>> cartSphCoeff;

// N    N is basisN if cartesian basis is used, sphBasisN if spherical basis is used

int N;

// cShellN                   number of contracted basis shells
// cShellLMax                maximum of cShellL
// cShellL                   contracted basis angular momentum
// cShellX                   x coordinate of contracted basis
// cShellY                   y coordinate of contracted basis
// cShellZ                   z coordinate of contracted basis
// cShellContraction         number of primitives in contraction
// cShellScale               scale of primitive function in contraction
// cShellPrimExp             exponent of primitive function in contraction
// cShellContractionCoeff    coefficient of primitive function in contraction

int cShellN = 0;
int cShellLMax;
std::vector<int> cShellL;
std::vector<real> cShellX;
std::vector<real> cShellY;
std::vector<real> cShellZ;
std::vector<int> cShellContraction;
std::vector<real> cShellScale;
std::vector<std::vector<real>> cShellPrimExp;
std::vector<std::vector<real>> cShellContractionCoeff;

// partitionAngularMomentum    given total angular momentum partitions it to cartesian component
// The following python script will generate these
// L = 
// ncart = (L + 1) * (L + 2) / 2;
// ang = []
// for i in range(L+1):
//     for j in range(i+1):
//         lx = L - i
//         ly = i - j
//         lz = j
//         ang.append([lx, ly, lz])

const std::vector<std::vector<std::vector<int>>> partitionAngularMomentum = {
    {{0, 0, 0}},
    {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}},
    {{2, 0, 0}, {1, 1, 0}, {1, 0, 1}, {0, 2, 0}, {0, 1, 1}, {0, 0, 2}},
    {{3, 0, 0}, {2, 1, 0}, {2, 0, 1}, {1, 2, 0}, {1, 1, 1}, {1, 0, 2}, {0, 3, 0}, {0, 2, 1}, {0, 1, 2}, {0, 0, 3}},
    {{4, 0, 0}, {3, 1, 0}, {3, 0, 1}, {2, 2, 0}, {2, 1, 1}, {2, 0, 2}, {1, 3, 0}, {1, 2, 1}, {1, 1, 2}, {1, 0, 3},
          {0, 4, 0}, {0, 3, 1}, {0, 2, 2}, {0, 1, 3}, {0, 0, 4}},
    {{5, 0, 0}, {4, 1, 0}, {4, 0, 1}, {3, 2, 0}, {3, 1, 1}, {3, 0, 2}, {2, 3, 0}, {2, 2, 1}, {2, 1, 2}, {2, 0, 3},
          {1, 4, 0}, {1, 3, 1}, {1, 2, 2}, {1, 1, 3}, {1, 0, 4}, {0, 5, 0}, {0, 4, 1}, {0, 3, 2}, {0, 2, 3}, {0, 1, 4},
          {0, 0, 5}},
    {{6, 0, 0}, {5, 1, 0}, {5, 0, 1}, {4, 2, 0}, {4, 1, 1}, {4, 0, 2}, {3, 3, 0}, {3, 2, 1}, {3, 1, 2}, {3, 0, 3},
          {2, 4, 0}, {2, 3, 1}, {2, 2, 2}, {2, 1, 3}, {2, 0, 4}, {1, 5, 0}, {1, 4, 1}, {1, 3, 2}, {1, 2, 3}, {1, 1, 4},
          {1, 0, 5}, {0, 6, 0}, {0, 5, 1}, {0, 4, 2}, {0, 3, 3}, {0, 2, 4}, {0, 1, 5}, {0, 0, 6}}
};

// map from cartesian to spherical functions
// literature:
//    Schlegel, H. B. & Frisch, M. J. Transformation between Cartesian and pure spherical harmonic Gaussians. Int. J. Quantum Chem. 54, 83–87 (1995).
//    https://theochem.github.io/horton/2.0.1/tech_ref_gaussian_basis.html
//    Coefficients generated from QM/SCRIPTS/Generate_Spherical_to_Cartesian.ipynb
// row representation
const std::vector<std::vector<std::vector<int>>> partitionSphContraction = {
    { // S
        {0}
    },
    { // P
        {2},
        {0},
        {1}
    },
    { // D
        {0, 3, 5},
        {2},
        {4},
        {0, 3},
        {1}
    },
    { // F
        {2, 7, 9},
        {0, 3, 5},
        {1, 6, 8},
        {2, 7},
        {4},
        {0, 3},
        {1, 6}
    },
    { // G
        {0, 3, 5, 10, 12, 14},
        {2, 7, 9},
        {4, 11, 13},
        {0, 5, 10, 12},
        {1, 6, 8},
        {2, 7},
        {4, 11},
        {0, 3, 10},
        {1, 6}
    },
    { // H
        {2, 7, 9, 16, 18, 20},
        {0, 3, 5, 10, 12, 14},
        {1, 6, 8, 15, 17, 19},
        {2, 9, 16, 18},
        {4, 11, 13},
        {0, 3, 5, 10, 12},
        {1, 6, 8, 15, 17},
        {2, 7, 16},
        {4, 11},
        {0, 3, 10},
        {1, 6, 15}
    },
    { // I
        {0, 3, 5, 10, 12, 14, 21, 23, 25, 27},
        {2, 7, 9, 16, 18, 20},
        {4, 11, 13, 22, 24, 26},
        {0, 3, 5, 10, 14, 21, 23, 25},
        {1, 6, 8, 15, 17, 19},
        {2, 7, 9, 16, 18},
        {4, 11, 13, 22, 24},
        {0, 3, 5, 10, 12, 21, 23},
        {1, 8, 15, 17},
        {2, 7, 16},
        {4, 11, 22},
        {0, 3, 10, 21},
        {1, 6, 15}
    }
};

// row representation
const std::vector<std::vector<std::vector<real>>> partitionSphCoeff = {
    { // S
        {1.}
    },
    { // P
        {1.},
        {1.},
        {1.}
    },
    { // D
        {-1./2., -1./2., 1.},
        {1.},
        {1.},
        {sqrt(3.)/2., -sqrt(3.)/2.},
        {1.}
    },
    { // F
        {-3.*sqrt(5.)/10., -3.*sqrt(5.)/10., 1.},
        {-sqrt(6.)/4., -sqrt(30.)/20., sqrt(30.)/5.},
        {-sqrt(30.)/20., -sqrt(6.)/4., sqrt(30.)/5.},
        {sqrt(3.)/2., -sqrt(3.)/2.},
        {1.},
        {sqrt(10.)/4., -3.*sqrt(2.)/4.},
        {3.*sqrt(2.)/4., -sqrt(10.)/4.}
    },
    { // G
        {3./8., 3.*sqrt(105.)/140., -3.*sqrt(105.)/35., 3./8., -3.*sqrt(105.)/35., 1.},
        {-3.*sqrt(70.)/28., -3.*sqrt(14.)/28., sqrt(70.)/7.},
        {-3.*sqrt(14.)/28., -3.*sqrt(70.)/28., sqrt(70.)/7.},
        {-sqrt(5.)/4., 3.*sqrt(21.)/14., sqrt(5.)/4., -3.*sqrt(21.)/14.},
        {-sqrt(35.)/14., -sqrt(35.)/14., 3.*sqrt(7.)/7.},
        {sqrt(10.)/4., -3.*sqrt(2.)/4.},
        {3.*sqrt(2.)/4., -sqrt(10.)/4.},
        {sqrt(35.)/8., -3.*sqrt(3.)/4., sqrt(35.)/8.},
        {sqrt(5.)/2., -sqrt(5.)/2.}
    },
    { // H
        {5./8., sqrt(105.)/28., -5.*sqrt(21.)/21., 5./8., -5.*sqrt(21.)/21., 1.},
        {sqrt(15.)/8., sqrt(35.)/28., -3.*sqrt(35.)/14., sqrt(15.)/24., -3.*sqrt(7.)/14., sqrt(15.)/3.},
        {sqrt(15.)/24., sqrt(35.)/28., -3.*sqrt(7.)/14., sqrt(15.)/8., -3.*sqrt(35.)/14., sqrt(15.)/3.},
        {-sqrt(105.)/12., sqrt(5.)/2., sqrt(105.)/12., -sqrt(5.)/2.},
        {-sqrt(15.)/6., -sqrt(15.)/6., sqrt(15.)/3.},
        {-sqrt(70.)/16., sqrt(30.)/24., sqrt(30.)/6., sqrt(70.)/16., -sqrt(6.)/2.},
        {-sqrt(70.)/16., -sqrt(30.)/24., sqrt(6.)/2., sqrt(70.)/16., -sqrt(30.)/6.},
        {sqrt(35.)/8., -3.*sqrt(3.)/4., sqrt(35.)/8.},
        {sqrt(5.)/2., -sqrt(5.)/2.},
        {3.*sqrt(14.)/16., -5.*sqrt(6.)/8., 5.*sqrt(14.)/16.},
        {5.*sqrt(14.)/16., -5.*sqrt(6.)/8., 3.*sqrt(14.)/16.}
    },
    { // I
        {-5./16., -5.*sqrt(33.)/176., 15.*sqrt(33.)/88., -5.*sqrt(33.)/176., 9.*sqrt(385.)/308., -5.*sqrt(33.)/22.,-5./16., 15.*sqrt(33.)/88., -5.*sqrt(33.)/22., 1.},
        {5.*sqrt(231.)/88., 5.*sqrt(11.)/44., -5.*sqrt(55.)/22., 5.*sqrt(231.)/264., -5.*sqrt(11.)/22., sqrt(231.)/11.},
        {5.*sqrt(231.)/264., 5.*sqrt(11.)/44., -5.*sqrt(11.)/22., 5.*sqrt(231.)/88., -5.*sqrt(55.)/22., sqrt(231.)/11.},
        {sqrt(210.)/32., sqrt(770.)/352., -sqrt(770.)/22., -sqrt(770.)/352., sqrt(770.)/22., -sqrt(210.)/32., sqrt(770.)/22., -sqrt(770.)/22.},
        {sqrt(2310.)/176., 5.*sqrt(22.)/88., -sqrt(110.)/11., sqrt(2310.)/176., -sqrt(110.)/11., sqrt(2310.)/33.},
        {-3.*sqrt(2310.)/176., 3.*sqrt(110.)/88., 5.*sqrt(22.)/22., 3.*sqrt(2310.)/176., -3.*sqrt(110.)/22.},
        {-3*sqrt(2310.)/176., -3.*sqrt(110.)/88., 3.*sqrt(110.)/22., 3.*sqrt(2310.)/176., -5.*sqrt(22.)/22.},
        {-3.*sqrt(7.)/16., 5.*sqrt(231.)/176., 5.*sqrt(231.)/88., 5.*sqrt(231.)/176., -9.*sqrt(55.)/44., -3.*sqrt(7.)/16., 5.*sqrt(231.)/88.},
        {-3.*sqrt(77.)/44., 5.*sqrt(33.)/22., 3.*sqrt(77.)/44., -5.*sqrt(33.)/22.},
        {3.*sqrt(14.)/16., -5.*sqrt(6.)/8., 5.*sqrt(14.)/16.},
        {5.*sqrt(14.)/16., -5.*sqrt(6.)/8., 3.*sqrt(14.)/16.},
        {sqrt(462.)/32., -15.*sqrt(14.)/32., 15.*sqrt(14.)/32., -sqrt(462.)/32.},
        {3.*sqrt(42.)/16., -5.*sqrt(10.)/8., 3.*sqrt(42.)/16.}
    }
};

// column representation
const std::vector<std::vector<std::vector<int>>> partitionCartSphContraction = {
    { // S
        {0}
    },
    { // P
        {1},
        {2},
        {0},
    },
    { // D
        {0,3,},
        {4,},
        {1,},
        {0,3,},
        {2,},
        {0,},
    },
    { // F
        {1,5,},
        {2,6,},
        {0,3,},
        {1,5,},
        {4,},
        {1,},
        {2,6,},
        {0,3,},
        {2,},
        {0,},
    },
    { // G
        {0,3,7,},
        {4,8,},
        {1,5,},
        {0,7,},
        {2,6,},
        {0,3,},
        {4,8,},
        {1,5,},
        {4,},
        {1,},
        {0,3,7,},
        {2,6,},
        {0,3,},
        {2,},
        {0,},
    },
    { // H
        {1,5,9,},
        {2,6,10,},
        {0,3,7,},
        {1,5,9,},
        {4,8,},
        {1,5,},
        {2,6,10,},
        {0,7,},
        {2,6,},
        {0,3,},
        {1,5,9,},
        {4,8,},
        {1,5,},
        {4,},
        {1,},
        {2,6,10,},
        {0,3,7,},
        {2,6,},
        {0,3,},
        {2,},
        {0,},
    },
    { // I
        {0,3,7,11,},
        {4,8,12,},
        {1,5,9,},
        {0,3,7,11,},
        {2,6,10,},
        {0,3,7,},
        {4,12,},
        {1,5,9,},
        {4,8,},
        {1,5,},
        {0,3,7,11,},
        {2,6,10,},
        {0,7,},
        {2,6,},
        {0,3,},
        {4,8,12,},
        {1,5,9,},
        {4,8,},
        {1,5,},
        {4,},
        {1,},
        {0,3,7,11,},
        {2,6,10,},
        {0,3,7,},
        {2,6,},
        {0,3,},
        {2,},
        {0,},
    }
};

// column representation
const std::vector<std::vector<std::vector<real>>> partitionCartSphCoeff = {
    { // S
        {1.}
    },
    { // P
        {1.},
        {1.},
        {1.},
    },
    { // D
        {-1./2.,sqrt(3.)/2.,},
        {1.,},
        {1.,},
        {-1./2.,-sqrt(3.)/2.,},
        {1.,},
        {1.,},
    },
    { // F
        {-sqrt(6.)/4.,sqrt(10.)/4.,},
        {-sqrt(30.)/20.,3.*sqrt(2.)/4.,},
        {-3.*sqrt(5.)/10.,sqrt(3.)/2.,},
        {-sqrt(30.)/20.,-3.*sqrt(2.)/4.,},
        {1.,},
        {sqrt(30.)/5.,},
        {-sqrt(6.)/4.,-sqrt(10.)/4.,},
        {-3.*sqrt(5.)/10.,-sqrt(3.)/2.,},
        {sqrt(30.)/5.,},
        {1.,},
    },
    { // G
        {3./8.,-sqrt(5.)/4.,sqrt(35.)/8.,},
        {-sqrt(35.)/14.,sqrt(5.)/2.,},
        {-3.*sqrt(70.)/28.,sqrt(10.)/4.,},
        {3.*sqrt(105.)/140.,-3.*sqrt(3.)/4.,},
        {-3.*sqrt(14.)/28.,3.*sqrt(2.)/4.,},
        {-3.*sqrt(105.)/35.,3.*sqrt(21.)/14.,},
        {-sqrt(35.)/14.,-sqrt(5.)/2.,},
        {-3.*sqrt(14.)/28.,-3.*sqrt(2.)/4.,},
        {3.*sqrt(7.)/7.,},
        {sqrt(70.)/7.,},
        {3./8.,sqrt(5.)/4.,sqrt(35.)/8.,},
        {-3.*sqrt(70.)/28.,-sqrt(10.)/4.,},
        {-3.*sqrt(105.)/35.,-3.*sqrt(21.)/14.,},
        {sqrt(70.)/7.,},
        {1.,},
    },
    { // H
        {sqrt(15.)/8.,-sqrt(70.)/16.,3.*sqrt(14.)/16.,},
        {sqrt(15.)/24.,-sqrt(70.)/16.,5.*sqrt(14.)/16.,},
        {5./8.,-sqrt(105.)/12.,sqrt(35.)/8.,},
        {sqrt(35.)/28.,sqrt(30.)/24.,-5.*sqrt(6.)/8.,},
        {-sqrt(15.)/6.,sqrt(5.)/2.,},
        {-3.*sqrt(35.)/14.,sqrt(30.)/6.,},
        {sqrt(35.)/28.,-sqrt(30.)/24.,-5.*sqrt(6.)/8.,},
        {sqrt(105.)/28.,-3.*sqrt(3.)/4.,},
        {-3.*sqrt(7.)/14.,sqrt(6.)/2.,},
        {-5.*sqrt(21.)/21.,sqrt(5.)/2.,},
        {sqrt(15.)/24.,sqrt(70.)/16.,5.*sqrt(14.)/16.,},
        {-sqrt(15.)/6.,-sqrt(5.)/2.,},
        {-3.*sqrt(7.)/14.,-sqrt(6.)/2.,},
        {sqrt(15.)/3.,},
        {sqrt(15.)/3.,},
        {sqrt(15.)/8.,sqrt(70.)/16.,3.*sqrt(14.)/16.,},
        {5./8.,sqrt(105.)/12.,sqrt(35.)/8.,},
        {-3.*sqrt(35.)/14.,-sqrt(30.)/6.,},
        {-5.*sqrt(21.)/21.,-sqrt(5.)/2.,},
        {sqrt(15.)/3.,},
        {1.,},
    },
    { // I
        {-5./16.,sqrt(210.)/32.,-3.*sqrt(7.)/16.,sqrt(462.)/32.,},
        {sqrt(2310.)/176.,-3.*sqrt(77.)/44.,3.*sqrt(42.)/16.,},
        {5.*sqrt(231.)/88.,-3.*sqrt(2310.)/176.,3.*sqrt(14.)/16.,},
        {-5.*sqrt(33.)/176.,sqrt(770.)/352.,5.*sqrt(231.)/176.,-15.*sqrt(14.)/32.,},
        {5.*sqrt(231.)/264.,-3.*sqrt(2310.)/176.,5.*sqrt(14.)/16.,},
        {15.*sqrt(33.)/88.,-sqrt(770.)/22.,5.*sqrt(231.)/88.,},
        {5.*sqrt(22.)/88.,-5.*sqrt(10.)/8.,},
        {5.*sqrt(11.)/44.,3.*sqrt(110.)/88.,-5.*sqrt(6.)/8.,},
        {-sqrt(110.)/11.,5.*sqrt(33.)/22.,},
        {-5.*sqrt(55.)/22.,5.*sqrt(22.)/22.,},
        {-5.*sqrt(33.)/176.,-sqrt(770.)/352.,5.*sqrt(231.)/176.,15.*sqrt(14.)/32.,},
        {5.*sqrt(11.)/44.,-3.*sqrt(110.)/88.,-5.*sqrt(6.)/8.,},
        {9.*sqrt(385.)/308.,-9.*sqrt(55.)/44.,},
        {-5.*sqrt(11.)/22.,3.*sqrt(110.)/22.,},
        {-5.*sqrt(33.)/22.,sqrt(770.)/22.,},
        {sqrt(2310.)/176.,3.*sqrt(77.)/44.,3.*sqrt(42.)/16.,},
        {5.*sqrt(231.)/264.,3.*sqrt(2310.)/176.,5.*sqrt(14.)/16.,},
        {-sqrt(110.)/11.,-5.*sqrt(33.)/22.,},
        {-5.*sqrt(11.)/22.,-3.*sqrt(110.)/22.,},
        {sqrt(2310.)/33.,},
        {sqrt(231.)/11.,},
        {-5./16.,-sqrt(210.)/32.,-3.*sqrt(7.)/16.,-sqrt(462.)/32.,},
        {5.*sqrt(231.)/88.,3.*sqrt(2310.)/176.,3.*sqrt(14.)/16.,},
        {15.*sqrt(33.)/88.,sqrt(770.)/22.,5.*sqrt(231.)/88.,},
        {-5.*sqrt(55.)/22.,-5.*sqrt(22.)/22.,},
        {-5.*sqrt(33.)/22.,-sqrt(770.)/22.,},
        {sqrt(231.)/11.,},
        {1.,},
    }
};

void normalizeContraction();
void buildSphContractionCoeff();
void buildCartSphContractionCoeff();


///////////////////////////////////////////////////////////////
//                                                           //
//  int sToL  --  return angular momentum for given orbital  //
//                                                           //
///////////////////////////////////////////////////////////////

int sToL(std::string orbital)
{
    if (orbital == "S")
        return 0;
    else if (orbital == "P")
        return 1;
    else if (orbital == "D")
        return 2;
    else if (orbital == "F")
        return 3;
    else if (orbital == "G")
        return 4;
    else if (orbital == "H")
        return 5;
    else if (orbital == "I")
        return 6;
    else
    {
        std::cerr << "Orbital must be one of S, P, D, F, G, H, or I." << std::endl;
        std::exit(1);
    }
}


/////////////////////////////////////
//                                 //
//  void kbasis  --  set up basis  //
//                                 //
/////////////////////////////////////

void kbasis()
{
    for (int i = 0; i < atoms::n; ++i)
    {
        std::string atomName = atoms::atom[i];
        real coordx = atoms::coordx[i];
        real coordy = atoms::coordy[i];
        real coordz = atoms::coordz[i];
        gbs::AoBasis &basis = gbs::AoBasisMap[atomName];
        int shellN = basis.getShellN();
        cShellN += shellN;
        std::vector<std::string>& shell = basis.getShell();
        std::vector<int>& contraction = basis.getContraction();
        std::vector<real>& scale = basis.getScale();
        std::vector<std::vector<real>>& primExp = basis.getPrimExp();
        std::vector<std::vector<real>>& contractionCoeff = basis.getContractionCoeff();
        for (int j = 0; j < shellN; ++j)
        {
            int l = sToL(shell[j]);
            basisN += lToN(l);
            sphBasisN += lToSphN(l);
            cShellL.push_back(l);
            cShellX.push_back(coordx);
            cShellY.push_back(coordy);
            cShellZ.push_back(coordz);
            cShellContraction.push_back(contraction[j]);
            cShellScale.push_back(scale[j]);
            cShellPrimExp.push_back(primExp[j]);
            cShellContractionCoeff.push_back(contractionCoeff[j]);
        }  
    }
    cShellLMax = *std::max_element(cShellL.begin(), cShellL.end());

    // set N to basisN or sphBasisN
    N = (gbs::basisType == gbs::BasisType::cartesian) ? basisN : sphBasisN;

    normalizeContraction();
    // // print to debug
    // std::cout << "Contraction normalization: " << std::endl;
    // for (int i = 0; i < basisN; ++i)
    // {
    //     printf("%20.16f", basisNorm[i]);
    // }
    // std::cout << std::endl;

    if (gbs::basisType == gbs::BasisType::spherical)
    {
        buildSphContractionCoeff();
        buildCartSphContractionCoeff();
    }
}


////////////////////////////////////////////////////////////
//                                                        //
//  void normalizeContraction  --  normalize contraction  //
//                                                        //
////////////////////////////////////////////////////////////

// Compute normalization Eq. 2.25 (Fermann, J. T. & Valeev, E. F. Fundamentals of Molecular Integrals Evaluation. (2020).)
// Also consult Cartesian Overlap documentation for derivation
void normalizeContraction()
{
    basisNorm.reserve(basisN);
    auto& partitionL = partitionAngularMomentum;
    // loop over shell
    for (int i = 0; i < cShellN; ++i)
    {
        int l = cShellL[i];
        auto& primExp = cShellPrimExp[i];
        auto& coeff = cShellContractionCoeff[i];
        int coeffN = coeff.size();
        for (auto& part: partitionL[l])
        {
            int lx = part[0];
            int ly = part[1];
            int lz = part[2];
            real s = 0.;
            for (int j = 0; j < coeffN; ++j)
            {
                real aj = primExp[j];
                real cj = coeff[j];
                for (int k = 0; k < coeffN; ++k)
                {
                    real ak = primExp[k];
                    real ck = coeff[k];
                    real ajtak = aj * ak;
                    real ajpak = aj + ak;
                    s += cj * ck * pow(2.*sqrt(ajtak)/ajpak, l+1.5);
                }
            }
            basisNorm.push_back(1/sqrt(s));
        }
    }
}


///////////////////////////////////////////////////////////////////////////////////
//                                                                               //
//  void buildSphContractionCoeff  --  build spherical contraction coefficients  //
//                                                                               //
///////////////////////////////////////////////////////////////////////////////////

void buildSphContractionCoeff()
{
    // reserve array
    sphContraction.reserve(sphBasisN);
    sphCoeff.reserve(sphBasisN);

    int counter = 0;

    // loop over shell
    for (int i = 0; i < cShellN; ++i)
    {
        int l = cShellL[i];
        auto& partSphCont = partitionSphContraction[l];
        auto& partSphCoeff = partitionSphCoeff[l];
        int sphN = lToSphN(l);
        for (int j = 0; j < sphN; ++j)
        {
            auto contraction = partSphCont[j];
            for (int k = 0; k < contraction.size(); ++k)
            {
                contraction[k] += counter;
            }
            sphContraction.push_back(contraction);
            sphCoeff.push_back(partSphCoeff[j]);
        }
        counter += lToN(l);
    }

    // // print to debug
    // for (int i = 0; i < sphBasisN; ++i)
    // {
    //     for (int j = 0; j < sphContraction[i].size(); ++j)
    //     {
    //         printf("%4i,", sphContraction[i][j]);
    //     }
    //     std::cout << std::endl;
    // }
    // std::cout << std::endl;
    // for (int i = 0; i < sphBasisN; ++i)
    // {
    //     for (int j = 0; j < sphCoeff[i].size(); ++j)
    //     {
    //         printf("%10.5f,", sphCoeff[i][j]);
    //     }
    //     std::cout << std::endl;
    // }
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                             //
//  void buildCartSphContractionCoeff  --  build spherical contraction coefficients for given cartesian basis  //
//                                                                                                             //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void buildCartSphContractionCoeff()
{
    // reserve array
    cartSphContraction.reserve(basisN);
    cartSphCoeff.reserve(basisN);

    int counter = 0;

    // loop over shell
    for (int i = 0; i < cShellN; ++i)
    {
        int l = cShellL[i];
        auto& partCont = partitionCartSphContraction[l];
        auto& partCoeff = partitionCartSphCoeff[l];
        int N = lToN(l);
        for (int j = 0; j < N; ++j)
        {
            auto contraction = partCont[j];
            for (int k = 0; k < contraction.size(); ++k)
            {
                contraction[k] += counter;
            }
            cartSphContraction.push_back(contraction);
            cartSphCoeff.push_back(partCoeff[j]);
        }
        counter += lToSphN(l);
    }

    // // print to debug
    // for (int i = 0; i < basisN; ++i)
    // {
    //     for (int j = 0; j < cartSphContraction[i].size(); ++j)
    //     {
    //         printf("%4i,", cartSphContraction[i][j]);
    //     }
    //     std::cout << std::endl;
    // }
    // std::cout << std::endl;
    // for (int i = 0; i < basisN; ++i)
    // {
    //     for (int j = 0; j < cartSphCoeff[i].size(); ++j)
    //     {
    //         printf("%10.5f,", cartSphCoeff[i][j]);
    //     }
    //     std::cout << std::endl;
    // }
}
}
