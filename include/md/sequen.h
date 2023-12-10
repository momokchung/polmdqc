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
//  sequen  --  sequence information for biopolymer  //
//                                                   //
///////////////////////////////////////////////////////

// nseq     total number of residues in biopolymer sequences
// nchain   number of separate biopolymer sequence chains
// ichain   first and last residue in each biopolymer chain
// seqtyp   residue type for each residue in the sequence
// seq      three-letter code for each residue in the sequence
// chnnam   one-letter identifier for each sequence chain
// chntyp   contents of each chain (GENERIC, PEPTIDE or NUCLEIC)

MDQC_EXTERN int nseq;
MDQC_EXTERN int nchain;
MDQC_EXTERN int ichain[maxres][2];
MDQC_EXTERN int seqtyp[maxres];
MDQC_EXTERN std::string chnnam[maxres];
MDQC_EXTERN std::string seq[maxres];
MDQC_EXTERN std::string chntyp[maxres];
}
