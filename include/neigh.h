/////////////////////////////////////////////////////////////
//                                                         //
//  neigh.h  --  pairwise neighbor list indices & storage  //
//                                                         //
/////////////////////////////////////////////////////////////

// maxvlst     maximum size of van der Waals pair neighbor lists
// maxelst     maximum size of electrostatic pair neighbor lists
// maxulst     maximum size of dipole preconditioner pair lists
// nvlst       number of sites in list for each vdw site
// vlst        site numbers in neighbor list of each vdw site
// nelst       number of sites in list for each electrostatic site
// elst        site numbers in list of each electrostatic site
// nulst       number of sites in list for each preconditioner site
// ulst        site numbers in list of each preconditioner site
// lbuffer     width of the neighbor list buffer region
// pbuffer     width of the preconditioner list buffer region
// lbuf2       square of half the neighbor list buffer width
// pbuf2       square of half the preconditioner list buffer width
// vbuf2       square of van der Waals cutoff plus the list buffer
// vbufx       square of van der Waals cutoff plus 2X list buffer
// dbuf2       square of dispersion cutoff plus the list buffer
// dbufx       square of dispersion cutoff plus 2X list buffer
// cbuf2       square of charge cutoff plus the list buffer
// cbufx       square of charge cutoff plus 2X list buffer
// mbuf2       square of multipole cutoff plus the list buffer
// mbufx       square of multipole cutoff plus 2X list buffer
// ubuf2       square of preconditioner cutoff plus the list buffer
// ubufx       square of preconditioner cutoff plus 2X list buffer
// xvold       x-coordinate at last vdw/dispersion list update
// yvold       y-coordinate at last vdw/dispersion list update
// zvold       z-coordinate at last vdw/dispersion list update
// xeold       x-coordinate at last electrostatic list update
// yeold       y-coordinate at last electrostatic list update
// zeold       z-coordinate at last electrostatic list update
// xuold       x-coordinate at last preconditioner list update
// yuold       y-coordinate at last preconditioner list update
// zuold       z-coordinate at last preconditioner list update
// dovlst      logical flag to rebuild vdw neighbor list
// dodlst      logical flag to rebuild dispersion neighbor list
// doclst      logical flag to rebuild charge neighbor list
// domlst      logical flag to rebuild multipole neighbor list
// doulst      logical flag to rebuild preconditioner neighbor list


#pragma once
#include "macro.h"
#include <string>
#include <vector>

QCMD_EXTERN int maxvlst;
QCMD_EXTERN int maxelst;
QCMD_EXTERN int maxulst;
QCMD_EXTERN std::vector<int> nvlst;
QCMD_EXTERN std::vector<std::vector<int>>vlst;
QCMD_EXTERN std::vector<int> nelst;
QCMD_EXTERN std::vector<std::vector<int>>elst;
QCMD_EXTERN std::vector<int> nulst;
QCMD_EXTERN std::vector<std::vector<int>>ulst;
QCMD_EXTERN double lbuffer,pbuffer;
QCMD_EXTERN double lbuf2,pbuf2;
QCMD_EXTERN double vbuf2,vbufx;
QCMD_EXTERN double dbuf2,dbufx;
QCMD_EXTERN double cbuf2,cbufx;
QCMD_EXTERN double mbuf2,mbufx;
QCMD_EXTERN double ubuf2,ubufx;
QCMD_EXTERN std::vector<double> xvold;
QCMD_EXTERN std::vector<double> yvold;
QCMD_EXTERN std::vector<double> zvold;
QCMD_EXTERN std::vector<double> xeold;
QCMD_EXTERN std::vector<double> yeold;
QCMD_EXTERN std::vector<double> zeold;
QCMD_EXTERN std::vector<double> xuold;
QCMD_EXTERN std::vector<double> yuold;
QCMD_EXTERN std::vector<double> zuold;
QCMD_EXTERN bool dovlst,dodlst;
QCMD_EXTERN bool doclst,domlst;
QCMD_EXTERN bool doulst;
 