// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include <string>

namespace polmdqc
{
/////////////////////////////////////////////////////////////
//                                                         //
//  cutoffSwitch  --  get switching function coefficients  //
//                                                         //
/////////////////////////////////////////////////////////////

enum class CutoffMode
{
    VdW,
    Repuls,
    Disp,
    Charge,
    ChgDpl,
    Dipole,
    Mpole,
    ChgTrn,
    Ewald,
    DEwald,
    USolv,
    GKV,
    GKSA,
};

void cutoffSwitch(CutoffMode mode);
}
