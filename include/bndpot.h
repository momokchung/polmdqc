//////////////////////////////////////////////////////////
//                                                      //
//  bndpot.h  --  bond stretch functional form details  //
//                                                      //
//////////////////////////////////////////////////////////

// cbnd      cubic coefficient in bond stretch potential
// qbnd      quartic coefficient in bond stretch potential
// bndunit   convert bond stretch energy to kcal/mole
// bndtyp    type of bond stretch potential energy function


#pragma once
#include "macro.h"
#include <string>

QCMD_EXTERN double cbnd;
QCMD_EXTERN double qbnd;
QCMD_EXTERN double bndunit;
QCMD_EXTERN std::string bndtyp;
