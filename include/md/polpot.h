// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"
#include <string>

namespace polmdqc
{
////////////////////////////////////////////////////////
//                                                    //
//  polpot  --  polarization functional form details  //
//                                                    //
////////////////////////////////////////////////////////

// politer      maximum number of induced dipole SCF iterations
// poleps       induced dipole convergence criterion (rms Debye/atom)
// p2scale      scale factor for 1-2 polarization energy interactions
// p3scale      scale factor for 1-3 polarization energy interactions
// p4scale      scale factor for 1-4 polarization energy interactions
// p5scale      scale factor for 1-5 polarization energy interactions
// p2iscale     scale factor for 1-2 intragroup polarization energy
// p3iscale     scale factor for 1-3 intragroup polarization energy
// p4iscale     scale factor for 1-4 intragroup polarization energy
// p5iscale     scale factor for 1-5 intragroup polarization energy
// d1scale      scale factor for intra-group direct induction
// d2scale      scale factor for 1-2 group direct induction
// d3scale      scale factor for 1-3 group direct induction
// d4scale      scale factor for 1-4 group direct induction
// u1scale      scale factor for intra-group mutual induction
// u2scale      scale factor for 1-2 group mutual induction
// u3scale      scale factor for 1-3 group mutual induction
// u4scale      scale factor for 1-4 group mutual induction
// w2scale      scale factor for 1-2 induced dipole interactions
// w3scale      scale factor for 1-3 induced dipole interactions
// w4scale      scale factor for 1-4 induced dipole interactions
// w5scale      scale factor for 1-5 induced dipole interactions
// uaccel       acceleration factor for induced dipole SCF iterations
// polprt       flag to print summary of induced dipole iterations
// dpequal      flag to set dscale values equal to pscale values
// use_thole    flag to use Thole damped polarization interactions
// use_tholed   flag to use alternate Thole for direct polarization
// use_expol    flag to use damped exchange polarization correction
// scrtyp       type of exchange polarization (S2U, S2 or G)
// poltyp       type of polarization (MUTUAL, DIRECT, OPT or TCG)

MDQC_EXTERN int politer;
MDQC_EXTERN real poleps;
MDQC_EXTERN real p2scale,p3scale;
MDQC_EXTERN real p4scale,p5scale;
MDQC_EXTERN real p2iscale,p3iscale;
MDQC_EXTERN real p4iscale,p5iscale;
MDQC_EXTERN real d1scale,d2scale;
MDQC_EXTERN real d3scale,d4scale;
MDQC_EXTERN real u1scale,u2scale;
MDQC_EXTERN real u3scale,u4scale;
MDQC_EXTERN real w2scale,w3scale;
MDQC_EXTERN real w4scale,w5scale;
MDQC_EXTERN real uaccel;
MDQC_EXTERN bool polprt;
MDQC_EXTERN bool dpequal;
MDQC_EXTERN bool use_thole;
MDQC_EXTERN bool use_tholed;
MDQC_EXTERN bool use_expol;
MDQC_EXTERN std::string scrtyp;
MDQC_EXTERN std::string poltyp;
}
