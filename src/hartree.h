///////////////////////////////////////////
//                                       //
//  hartree.h  --  compute Hartree Fock  //
//                                       //
///////////////////////////////////////////


#pragma once
#include "mathUtils.h"

namespace hartree
{
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

void eigenS(std::vector<std::vector<real>>& S);
void guess();
void rhf();
void scf();
}
