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

MDQC_EXTERN double angunit,stbnunit;
MDQC_EXTERN double aaunit,opbunit;
MDQC_EXTERN double opdunit;
MDQC_EXTERN double cang,qang;
MDQC_EXTERN double pang,sang;
MDQC_EXTERN double copb,qopb;
MDQC_EXTERN double popb,sopb;
MDQC_EXTERN double copd,qopd;
MDQC_EXTERN double popd,sopd;
MDQC_EXTERN std::string opbtyp;
MDQC_EXTERN std::vector<std::string> angtyp;
}