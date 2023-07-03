///////////////////////////////
//                           //
//  kbasis.h  --  basis set  //
//                           //
///////////////////////////////


#pragma once
#include "init.h"
#include <string>
#include <vector>

namespace basis
{
extern int basisN;
// extern std::vector<int> basisLx;
// extern std::vector<int> basisLy;
// extern std::vector<int> basisLz;
// extern std::vector<int> basisL;
// extern std::vector< std::vector<real> > basisCoeff;
// extern std::vector< std::vector<real> > Exp;
// extern std::vector<real> basisX;
// extern std::vector<real> basisY;
// extern std::vector<real> basisZ;
extern std::vector<real> basisNorm;

extern int N;

extern int sphBasisN;
extern std::vector< std::vector<int> > sphContraction;
extern std::vector< std::vector<real> > sphCoeff;
extern std::vector< std::vector<int> > cartSphContraction;
extern std::vector< std::vector<real> > cartSphCoeff;

extern int primN;
extern std::vector<real> primNorm;
extern std::vector<int> primToBasis;

extern int cShellN;
extern int cShellLMax;
extern std::vector<int> cShellL;
extern std::vector<real> cShellX;
extern std::vector<real> cShellY;
extern std::vector<real> cShellZ;
extern std::vector<int> cShellContraction;
extern std::vector<real> cShellScale;
extern std::vector< std::vector<real> > cShellPrimExp;
extern std::vector< std::vector<real> > cShellContractionCoeff;

extern const std::vector< std::vector< std::vector<int> > > partitionAngularMomentum;

void kbasis();


//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//  inline int lToN  --  return number of basis for given angular momentum  //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

inline int lToN(int l)
{
    return (l + 1) * (l + 2) / 2;
}


///////////////////////////////////////////////////////////////////////////////////////////
//                                                                                       //
//  inline int lToSphN  --  return number of spherical basis for given angular momentum  //
//                                                                                       //
///////////////////////////////////////////////////////////////////////////////////////////

inline int lToSphN(int l)
{
    return 2 * l + 1;
}
}
