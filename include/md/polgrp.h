// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"
#include <vector>

namespace polmdqc
{
/////////////////////////////////////////////////////////
//                                                     //
//  polgrp  --  polarization group connectivity lists  //
//                                                     //
/////////////////////////////////////////////////////////

// maxp11   maximum number of atoms in a polarization group
// maxp12   maximum number of atoms in groups 1-2 to an atom
// maxp13   maximum number of atoms in groups 1-3 to an atom
// maxp14   maximum number of atoms in groups 1-4 to an atom
// np11     number of atoms in polarization group of each atom
// np12     number of atoms in groups 1-2 to each atom
// np13     number of atoms in groups 1-3 to each atom
// np14     number of atoms in groups 1-4 to each atom
// ip11     atom numbers of atoms in same group as each atom
// ip12     atom numbers of atoms in groups 1-2 to each atom
// ip13     atom numbers of atoms in groups 1-3 to each atom
// ip14     atom numbers of atoms in groups 1-4 to each atom

constexpr int maxp11 = 200;
constexpr int maxp12 = 200;
constexpr int maxp13 = 200;
constexpr int maxp14 = 200;
MDQC_EXTERN std::vector<int> np11;
MDQC_EXTERN std::vector<int> np12;
MDQC_EXTERN std::vector<int> np13;
MDQC_EXTERN std::vector<int> np14;
MDQC_EXTERN std::vector<std::vector<int>> ip11;
MDQC_EXTERN std::vector<std::vector<int>> ip12;
MDQC_EXTERN std::vector<std::vector<int>> ip13;
MDQC_EXTERN std::vector<std::vector<int>> ip14;
}
