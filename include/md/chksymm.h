// Author: Moses KJ Chung
// Year:   2024

#pragma once
#include "precision.h"
#include <string>

namespace polmdqc
{
////////////////////////////////////////////////////////////////
//                                                            //
//  subroutine chksymm  --  test for 1D, 2D & other symmetry  //
//                                                            //
////////////////////////////////////////////////////////////////

enum class SymTyp
{
    None,
    Single,
    Linear,
    Planar,
    Mirror,
    Center,
};

void chksymm(int n, real* mass, real* xref, real* yref, real* zref, SymTyp& symtyp);
}
