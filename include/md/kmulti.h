// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"
#include <string>
#include <vector>

namespace polmdqc
{
//////////////////////////////////////////////////////////
//                                                      //
//  kmulti  --  atomic multipole forcefield parameters  //
//                                                      //
//////////////////////////////////////////////////////////

// LocalFrame   local frame enum class
// maxnmp       maximum number of atomic multipole parameter entries
// multip       atomic monopole, dipole and quadrupole values
// mpaxis       type of local axis definition for atomic multipoles
// kmp          string of atom types for atomic multipoles
// retLFRM      returns local frame

enum class LocalFrame
{
    None,
    ZthenX,
    ZOnly,
    Bisector,
    ZBisect,
    ThreeFold,
};

MDQC_EXTERN int maxnmp;
MDQC_EXTERN std::vector<std::vector<real>> multip;
MDQC_EXTERN std::vector<LocalFrame> mpaxis;
MDQC_EXTERN std::vector<std::string> kmp;

inline LocalFrame retLFRM(std::string axt)
{
    LocalFrame lfrm = LocalFrame::None;
    if (axt == "Z-then-X") lfrm = LocalFrame::ZthenX;
    else if (axt == "Z-Only") lfrm = LocalFrame::ZOnly;
    else if (axt == "Bisector") lfrm = LocalFrame::Bisector;
    else if (axt == "Z-Bisect") lfrm = LocalFrame::ZBisect;
    else if (axt == "3-Fold") lfrm = LocalFrame::ThreeFold;
    return lfrm;
}
}
