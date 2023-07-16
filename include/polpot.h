//////////////////////////////////////////////////////////
//                                                      //
//  polpot.h  --  polarization functional form details  //
//                                                      //
//////////////////////////////////////////////////////////

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


#pragma once
#include "macro.h"
#include <string>

QCMD_EXTERN int politer;
QCMD_EXTERN double poleps;
QCMD_EXTERN double p2scale,p3scale;
QCMD_EXTERN double p4scale,p5scale;
QCMD_EXTERN double p2iscale,p3iscale;
QCMD_EXTERN double p4iscale,p5iscale;
QCMD_EXTERN double d1scale,d2scale;
QCMD_EXTERN double d3scale,d4scale;
QCMD_EXTERN double u1scale,u2scale;
QCMD_EXTERN double u3scale,u4scale;
QCMD_EXTERN double w2scale,w3scale;
QCMD_EXTERN double w4scale,w5scale;
QCMD_EXTERN double uaccel;
QCMD_EXTERN bool polprt;
QCMD_EXTERN bool dpequal;
QCMD_EXTERN bool use_thole;
QCMD_EXTERN bool use_tholed;
QCMD_EXTERN bool use_expol;
QCMD_EXTERN std::string scrtyp;
QCMD_EXTERN std::string poltyp;
