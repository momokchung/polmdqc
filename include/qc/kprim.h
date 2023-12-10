// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "init.h"
#include <vector>

namespace prim
{
//////////////////////////////////////
//                                  //
//  kprim  --  primitive gaussians  //
//                                  //
//////////////////////////////////////

extern int primN;
// extern std::vector<int> primLx;
// extern std::vector<int> primLy;
// extern std::vector<int> primLz;
// extern std::vector<real> primExp;
// extern std::vector<real> primX;
// extern std::vector<real> primY;
// extern std::vector<real> primZ;
// extern std::vector<real> primScale;
extern std::vector<real> primNorm;
extern std::vector<int> primToBasis;
extern std::vector<int> primToShell;

extern int primShellN;
extern std::vector<int> primShellL;
extern std::vector<int> primShellIndex;
extern std::vector<real> primShellX;
extern std::vector<real> primShellY;
extern std::vector<real> primShellZ;
extern std::vector<real> primShellExp;
extern std::vector<real> primShellCoeff;
extern std::vector<real> primShellScale;

void kprim();

}
