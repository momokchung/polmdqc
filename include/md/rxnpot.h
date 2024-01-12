// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"

namespace polmdqc
{
//////////////////////////////////////////////////////////
//                                                      //
//  rxnpot  --  reaction field functional form details  //
//                                                      //
//////////////////////////////////////////////////////////

// rfsize    radius of reaction field sphere centered at origin
// rfbulkd   bulk dielectric constant of reaction field continuum
// rfterms   number of terms to use in reaction field summation

MDQC_EXTERN int rfterms;
MDQC_EXTERN real rfsize;
MDQC_EXTERN real rfbulkd;
}
