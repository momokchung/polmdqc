/////////////////////////////////////////////////////////
//                                                     //
//  sequen.h  --  sequence information for biopolymer  //
//                                                     //
/////////////////////////////////////////////////////////

// nseq     total number of residues in biopolymer sequences
// nchain   number of separate biopolymer sequence chains
// ichain   first and last residue in each biopolymer chain
// seqtyp   residue type for each residue in the sequence
// seq      three-letter code for each residue in the sequence
// chnnam   one-letter identifier for each sequence chain
// chntyp   contents of each chain (GENERIC, PEPTIDE or NUCLEIC)

#pragma once
#include "macro.h"
#include "sizes.h"
#include <string>

QCMD_EXTERN int nseq;
QCMD_EXTERN int nchain;
QCMD_EXTERN int ichain[maxres][2];
QCMD_EXTERN int seqtyp[maxres];
QCMD_EXTERN std::string chnnam[maxres];
QCMD_EXTERN std::string seq[maxres];
QCMD_EXTERN std::string chntyp[maxres];
