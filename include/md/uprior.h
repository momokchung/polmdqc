// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"
#include <string>
#include <vector>

namespace polmdqc
{
//////////////////////////////////////////////////////
//                                                  //
//  uprior  --  previous values of induced dipoles  //
//                                                  //
//////////////////////////////////////////////////////

// maxpred   maximum number of predictor induced dipoles to save
// nualt     number of sets of prior induced dipoles in storage
// maxualt   number of sets of induced dipoles needed for predictor
// gear      coefficients for Gear predictor binomial method
// aspc      coefficients for always stable predictor-corrector
// bpred     coefficients for induced dipole predictor polynomial
// bpredp    coefficients for predictor polynomial in energy field
// bpreds    coefficients for predictor for PB/GK solvation
// bpredps   coefficients for predictor in PB/GK energy field
// udalt     prior values for induced dipoles at each site
// upalt     prior values for induced dipoles in energy field
// usalt     prior values for induced dipoles for PB/GK solvation
// upsalt    prior values for induced dipoles in PB/GK energy field
// use_pred  flag to control use of induced dipole prediction
// polpred   type of predictor polynomial (ASPC, GEAR or LSQR)

constexpr int maxpred = 17;
MDQC_EXTERN int nualt;
MDQC_EXTERN int maxualt;
MDQC_EXTERN real gear[maxpred];
MDQC_EXTERN real aspc[maxpred];
MDQC_EXTERN real bpred[maxpred];
MDQC_EXTERN real bpredp[maxpred];
MDQC_EXTERN real bpreds[maxpred];
MDQC_EXTERN real bpredps[maxpred];
MDQC_EXTERN std::vector<std::vector<std::vector<real>>> udalt;
MDQC_EXTERN std::vector<std::vector<std::vector<real>>> upalt;
MDQC_EXTERN std::vector<std::vector<std::vector<real>>> usalt;
MDQC_EXTERN std::vector<std::vector<std::vector<real>>> upsalt;
MDQC_EXTERN bool use_pred;
MDQC_EXTERN std::string polpred;
}
