// Author: Moses KJ Chung
// Year:   2024

#pragma once
#include "macro.h"

namespace polmdqc
{
////////////////////////////////////////////////////
//                                                //
//  methqm  --  theory/method for QM calculation  //
//                                                //
////////////////////////////////////////////////////

// QMMethodType   enum class for QM method
// qmmethod       QM method to use for calculation

enum class QMMethodType
{
    HF,
    MP2,
    CCSD,
    CCSDt,
    CCSDT,
    SAPT0,
    SAPT2,
    SAPT2p,
    ALMO,
};
MDQC_EXTERN QMMethodType qmmethod;
}
