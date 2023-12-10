// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "mathUtils.h"

namespace hartree
{
/////////////////////////////////////////
//                                     //
//  hartree  --  compute Hartree Fock  //
//                                     //
/////////////////////////////////////////

enum class SCFType
{
    incore,
    conventional,
    direct,
    densityFit,
};

enum class SCFConvergence
{
    core,
    diis
};

// set s_tolerance from xyz/key file in the future
extern real s_tolerance;

extern std::vector<real> D;
extern std::vector<real> G;

extern std::vector<std::vector<real>> S;
extern std::vector<std::vector<real>> KE;
extern std::vector<std::vector<real>> NE;

void eigenS(std::vector<std::vector<real>>& S);
void guess();
void rhf();
void scf();
}
