// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"
#include <string>

namespace polmdqc
{
//////////////////////////////////////////////////////////
//                                                      //
//  resdue  --  amino acid & nucleotide residue names   //
//                                                      //
//////////////////////////////////////////////////////////

// maxamino  maximum number of amino acid residue types
// maxnuc    maximum number of nucleic acid residue types
// ntyp      biotypes for mid-chain peptide backbone N atoms
// catyp     biotypes for mid-chain peptide backbone CA atoms
// ctyp      biotypes for mid-chain peptide backbone C atoms
// hntyp     biotypes for mid-chain peptide backbone HN atoms
// otyp      biotypes for mid-chain peptide backbone O atoms
// hatyp     biotypes for mid-chain peptide backbone HA atoms
// cbtyp     biotypes for mid-chain peptide backbone CB atoms
// nntyp     biotypes for N-terminal peptide backbone N atoms
// cantyp    biotypes for N-terminal peptide backbone CA atoms
// cntyp     biotypes for N-terminal peptide backbone C atoms
// hnntyp    biotypes for N-terminal peptide backbone HN atoms
// ontyp     biotypes for N-terminal peptide backbone O atoms
// hantyp    biotypes for N-terminal peptide backbone HA atoms
// nctyp     biotypes for C-terminal peptide backbone N atoms
// cactyp    biotypes for C-terminal peptide backbone CA atoms
// cctyp     biotypes for C-terminal peptide backbone C atoms
// hnctyp    biotypes for C-terminal peptide backbone HN atoms
// octyp     biotypes for C-terminal peptide backbone O atoms
// hactyp    biotypes for C-terminal peptide backbone HA atoms
// o5typ     biotypes for nucleotide backbone and sugar O5' atoms
// c5typ     biotypes for nucleotide backbone and sugar C5' atoms
// h51typ    biotypes for nucleotide backbone and sugar H5' atoms
// h52typ    biotypes for nucleotide backbone and sugar H5'' atoms
// c4typ     biotypes for nucleotide backbone and sugar C4' atoms
// h4typ     biotypes for nucleotide backbone and sugar H4' atoms
// o4typ     biotypes for nucleotide backbone and sugar O4' atoms
// c1typ     biotypes for nucleotide backbone and sugar C1' atoms
// h1typ     biotypes for nucleotide backbone and sugar H1' atoms
// c3typ     biotypes for nucleotide backbone and sugar C3' atoms
// h3typ     biotypes for nucleotide backbone and sugar H3' atoms
// c2typ     biotypes for nucleotide backbone and sugar C2' atoms
// h21typ    biotypes for nucleotide backbone and sugar H2' atoms
// o2typ     biotypes for nucleotide backbone and sugar O2' atoms
// h22typ    biotypes for nucleotide backbone and sugar H2'' atoms
// o3typ     biotypes for nucleotide backbone and sugar O3' atoms
// ptyp      biotypes for nucleotide backbone and sugar P atoms
// optyp     biotypes for nucleotide backbone and sugar OP atoms
// h5ttyp    biotypes for nucleotide backbone and sugar H5T atoms
// h3ttyp    biotypes for nucleotide backbone and sugar H3T atoms
// amino     three-letter abbreviations for amino acids types
// nuclz     three-letter abbreviations for nucleic acids types
// amino1    one-letter abbreviations for amino acids types
// nuclz1    one-letter abbreviations for nucleic acids types

constexpr int maxamino=38;
constexpr int maxnuc=12;
MDQC_EXTERN int ntyp[maxamino];
MDQC_EXTERN int catyp[maxamino];
MDQC_EXTERN int ctyp[maxamino];
MDQC_EXTERN int hntyp[maxamino];
MDQC_EXTERN int otyp[maxamino];
MDQC_EXTERN int hatyp[maxamino];
MDQC_EXTERN int cbtyp[maxamino];
MDQC_EXTERN int nntyp[maxamino];
MDQC_EXTERN int cantyp[maxamino];
MDQC_EXTERN int cntyp[maxamino];
MDQC_EXTERN int hnntyp[maxamino];
MDQC_EXTERN int ontyp[maxamino];
MDQC_EXTERN int hantyp[maxamino];
MDQC_EXTERN int nctyp[maxamino];
MDQC_EXTERN int cactyp[maxamino];
MDQC_EXTERN int cctyp[maxamino];
MDQC_EXTERN int hnctyp[maxamino];
MDQC_EXTERN int octyp[maxamino];
MDQC_EXTERN int hactyp[maxamino];
MDQC_EXTERN int o5typ[maxnuc];
MDQC_EXTERN int c5typ[maxnuc];
MDQC_EXTERN int h51typ[maxnuc];
MDQC_EXTERN int h52typ[maxnuc];
MDQC_EXTERN int c4typ[maxnuc];
MDQC_EXTERN int h4typ[maxnuc];
MDQC_EXTERN int o4typ[maxnuc];
MDQC_EXTERN int c1typ[maxnuc];
MDQC_EXTERN int h1typ[maxnuc];
MDQC_EXTERN int c3typ[maxnuc];
MDQC_EXTERN int h3typ[maxnuc];
MDQC_EXTERN int c2typ[maxnuc];
MDQC_EXTERN int h21typ[maxnuc];
MDQC_EXTERN int o2typ[maxnuc];
MDQC_EXTERN int h22typ[maxnuc];
MDQC_EXTERN int o3typ[maxnuc];
MDQC_EXTERN int ptyp[maxnuc];
MDQC_EXTERN int optyp[maxnuc];
MDQC_EXTERN int h5ttyp[maxnuc];
MDQC_EXTERN int h3ttyp[maxnuc];
MDQC_EXTERN std::string amino1[maxamino];
MDQC_EXTERN std::string nuclz1[maxnuc];
MDQC_EXTERN std::string amino[maxamino];
MDQC_EXTERN std::string nuclz[maxnuc];
}
