// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"
#include <string>
#include <vector>

namespace polmdqc
{
///////////////////////////////////////////////////////////
//                                                       //
//  neigh  --  pairwise neighbor list indices & storage  //
//                                                       //
///////////////////////////////////////////////////////////

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

MDQC_EXTERN int maxvlst;
MDQC_EXTERN int maxelst;
MDQC_EXTERN int maxulst;
MDQC_EXTERN std::vector<int> nvlst;
MDQC_EXTERN std::vector<std::vector<int>>vlst;
MDQC_EXTERN std::vector<int> nelst;
MDQC_EXTERN std::vector<std::vector<int>>elst;
MDQC_EXTERN std::vector<int> nulst;
MDQC_EXTERN std::vector<std::vector<int>>ulst;
MDQC_EXTERN real lbuffer,pbuffer;
MDQC_EXTERN real lbuf2,pbuf2;
MDQC_EXTERN real vbuf2,vbufx;
MDQC_EXTERN real dbuf2,dbufx;
MDQC_EXTERN real cbuf2,cbufx;
MDQC_EXTERN real mbuf2,mbufx;
MDQC_EXTERN real ubuf2,ubufx;
MDQC_EXTERN std::vector<real> xvold;
MDQC_EXTERN std::vector<real> yvold;
MDQC_EXTERN std::vector<real> zvold;
MDQC_EXTERN std::vector<real> xeold;
MDQC_EXTERN std::vector<real> yeold;
MDQC_EXTERN std::vector<real> zeold;
MDQC_EXTERN std::vector<real> xuold;
MDQC_EXTERN std::vector<real> yuold;
MDQC_EXTERN std::vector<real> zuold;
MDQC_EXTERN bool dovlst,dodlst;
MDQC_EXTERN bool doclst,domlst;
MDQC_EXTERN bool doulst;
}
