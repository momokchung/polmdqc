// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"
#include "sizes.h"
#include <string>

namespace polmdqc
{
///////////////////////////////////////////////////////
//                                                   //
//  pdb  --  Protein Data Bank structure definition  //
//                                                   //
///////////////////////////////////////////////////////

// npdb      number of atoms stored in Protein Data Bank format
// nres      number of residues stored in Protein Data Bank format
// resnum    number of the residue to which each atom belongs
// resatm    number of first and last atom in each residue
// npdb12    number of atoms directly bonded to each CONECT atom
// ipdb12    atom numbers of atoms connected to each CONECT atom
// pdblist   list of the Protein Data Bank atom number of each atom
// xpdb      x-coordinate of each atom stored in PDB format
// ypdb      y-coordinate of each atom stored in PDB format
// zpdb      z-coordinate of each atom stored in PDB format
// altsym    string with PDB alternate locations to be included
// pdbres    Protein Data Bank residue name assigned to each atom
// pdbsym    Protein Data Bank atomic symbol assigned to each atom
// pdbatm    Protein Data Bank atom name assigned to each atom
// pdbtyp    Protein Data Bank record type assigned to each atom
// chnsym    string with PDB chain identifiers to be included
// instyp    string with PDB insertion records to be included

MDQC_EXTERN int npdb,nres;
MDQC_EXTERN MDQCArray<int> resnum;
MDQC_EXTERN MDQCArray2D<int,2> resatm;
MDQC_EXTERN MDQCArray<int> npdb12;
MDQC_EXTERN MDQCArray2D<int,maxval> ipdb12;
MDQC_EXTERN MDQCArray<int> pdblist;
MDQC_EXTERN MDQCArray<real> xpdb;
MDQC_EXTERN MDQCArray<real> ypdb;
MDQC_EXTERN MDQCArray<real> zpdb;
MDQC_EXTERN std::string altsym;
MDQC_EXTERN MDQCArray<std::string> pdbres;
MDQC_EXTERN MDQCArray<std::string> pdbsym;
MDQC_EXTERN MDQCArray<std::string> pdbatm;
MDQC_EXTERN MDQCArray<std::string> pdbtyp;
MDQC_EXTERN std::string chnsym;
MDQC_EXTERN std::string instyp;
}
