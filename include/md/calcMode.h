// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"

namespace polmdqc
{
///////////////////////////////////////////
//                                       //
//  calcMode  --  calculation mode type  //
//                                       //
///////////////////////////////////////////

// c     CalcMode   calculation modes
// c     CalcFlag   calculation flags

enum class CalcMode
{
    None,
    Energy,
    Analysis,
    Gradient,
    Virial,
    Hessian,
};

struct CalcFlag {
    bool do_energy;
    bool do_analysis;
    bool do_gradient;
    bool do_virial;
};

template <CalcMode mode>
constexpr CalcFlag getCalculationFlags() {
    switch (mode) {
        case CalcMode::None:
            return {false, false, false, false};
        case CalcMode::Energy:
            return {true, false, false, false};
        case CalcMode::Analysis:
            return {true, true, false, false};
        case CalcMode::Gradient:
            return {true, false, true, false};
        case CalcMode::Virial:
            return {true, false, true, true};
    }
}
}
