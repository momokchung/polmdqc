/////////////////////////////////////////////////////////
//                                                     //
//  pdb.h  --  Protein Data Bank structure definition  //
//                                                     //
/////////////////////////////////////////////////////////

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


#pragma once
#include "macro.h"
#include <string>
#include <vector>

QCMD_EXTERN int npdb,nres;
QCMD_EXTERN std::vector<int> resnum;
QCMD_EXTERN std::vector<std::vector<int>> resatm;
QCMD_EXTERN std::vector<int> npdb12;
QCMD_EXTERN std::vector<std::vector<int>> ipdb12;
QCMD_EXTERN std::vector<int> pdblist;
QCMD_EXTERN std::vector<double> xpdb;
QCMD_EXTERN std::vector<double> ypdb;
QCMD_EXTERN std::vector<double> zpdb;
QCMD_EXTERN std::string altsym;
QCMD_EXTERN std::vector<std::string> pdbres;
QCMD_EXTERN std::vector<std::string> pdbsym;
QCMD_EXTERN std::vector<std::string> pdbatm;
QCMD_EXTERN std::vector<std::string> pdbtyp;
QCMD_EXTERN std::string chnsym;
QCMD_EXTERN std::string instyp;
