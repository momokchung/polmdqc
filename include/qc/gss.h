// Author: Moses KJ Chung
// Year:   2024

#pragma once
#include "macro.h"

namespace polmdqc
{
/////////////////////////////////////////////
//                                         //
//  gss  --  guess parameters for quantum  //
//                                         //
/////////////////////////////////////////////

// DenGuess   enum class for initial guess of the density matrix
// bssguess   guess basis with small basis SCF projected into target basis
// denguess   initial guess of the density matrix


enum class DenGuess
{
    Core,
    SAD,
    SAP,
};
MDQC_EXTERN bool bssguess;
MDQC_EXTERN DenGuess denguess;
}
