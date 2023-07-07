///////////////////////////////////////////////////////
//                                                   //
//  sizes.h  --  parameters to set array dimensions  //
//                                                   //
///////////////////////////////////////////////////////

// "sizes" sets values for array dimensions used throughout
// the software; these parameters fix the size of the largest
// systems that can be handled

// maxatm          atoms in the molecular system
// maxtyp          force field atom type definitions
// maxclass        force field atom class definitions
// maxval          atoms directly bonded to an atom
// maxref          stored reference molecular systems
// maxgrp          user-defined groups of atoms
// maxres          residues in all macromolecules
// maxbio          biopolymer atom type definitions


#pragma once
#include <string>

constexpr int maxatm=1000000;
constexpr int maxtyp=5000;
constexpr int maxclass=1000;
constexpr int maxval=8;
constexpr int maxref=30;
constexpr int maxgrp=1000;
constexpr int maxres=10000;
constexpr int maxbio=10000;
