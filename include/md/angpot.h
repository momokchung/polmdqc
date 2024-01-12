// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"
#include <string>
#include <vector>

namespace polmdqc
{
//////////////////////////////////////////////////////
//                                                  //
//  angpot  --  angle bend functional form details  //
//                                                  //
//////////////////////////////////////////////////////

// angunit    convert angle bending energy to kcal/mole
// stbnunit   convert stretch-bend energy to kcal/mole
// aaunit     convert angle-angle energy to kcal/mole
// opbunit    convert out-of-plane bend energy to kcal/mole
// opdunit    convert out-of-plane distance energy to kcal/mole
// cang       cubic coefficient in angle bending potential
// qang       quartic coefficient in angle bending potential
// pang       quintic coefficient in angle bending potential
// sang       sextic coefficient in angle bending potential
// copb       cubic coefficient in out-of-plane bend potential
// qopb       quartic coefficient in out-of-plane bend potential
// popb       quintic coefficient in out-of-plane bend potential
// sopb       sextic coefficient in out-of-plane bend potential
// copd       cubic coefficient in out-of-plane distance potential
// qopd       quartic coefficient in out-of-plane distance potential
// popd       quintic coefficient in out-of-plane distance potential
// sopd       sextic coefficient in out-of-plane distance potential
// opbtyp     type of out-of-plane bend potential energy function
// angtyp     type of angle bending function for each bond angle

MDQC_EXTERN real angunit,stbnunit;
MDQC_EXTERN real aaunit,opbunit;
MDQC_EXTERN real opdunit;
MDQC_EXTERN real cang,qang;
MDQC_EXTERN real pang,sang;
MDQC_EXTERN real copb,qopb;
MDQC_EXTERN real popb,sopb;
MDQC_EXTERN real copd,qopd;
MDQC_EXTERN real popd,sopd;
MDQC_EXTERN std::string opbtyp;
MDQC_EXTERN std::vector<std::string> angtyp;
}
