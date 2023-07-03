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

const int maxatm=1000000;
const int maxtyp=5000;
const int maxclass=1000;
const int maxval=8;
const int maxref=30;
const int maxgrp=1000;
const int maxres=10000;
const int maxbio=10000;
