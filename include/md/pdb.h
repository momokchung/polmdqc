// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"
#include <string>
#include <vector>

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
MDQC_EXTERN std::vector<int> resnum;
MDQC_EXTERN std::vector<std::vector<int>> resatm;
MDQC_EXTERN std::vector<int> npdb12;
MDQC_EXTERN std::vector<std::vector<int>> ipdb12;
MDQC_EXTERN std::vector<int> pdblist;
MDQC_EXTERN std::vector<real> xpdb;
MDQC_EXTERN std::vector<real> ypdb;
MDQC_EXTERN std::vector<real> zpdb;
MDQC_EXTERN std::string altsym;
MDQC_EXTERN std::vector<std::string> pdbres;
MDQC_EXTERN std::vector<std::string> pdbsym;
MDQC_EXTERN std::vector<std::string> pdbatm;
MDQC_EXTERN std::vector<std::string> pdbtyp;
MDQC_EXTERN std::string chnsym;
MDQC_EXTERN std::string instyp;
}
