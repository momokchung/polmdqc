////////////////////////////////////////////////////////////
//                                                        //
//  rxnpot.h  --  reaction field functional form details  //
//                                                        //
////////////////////////////////////////////////////////////

// rfsize    radius of reaction field sphere centered at origin
// rfbulkd   bulk dielectric constant of reaction field continuum
// rfterms   number of terms to use in reaction field summation


#pragma once
#include "macro.h"

QCMD_EXTERN int rfterms;
QCMD_EXTERN double rfsize;
QCMD_EXTERN double rfbulkd;
