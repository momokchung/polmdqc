////////////////////////////////////////////////////////
//                                                    //
//  uprior.h  --  previous values of induced dipoles  //
//                                                    //
////////////////////////////////////////////////////////

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


#pragma once
#include "macro.h"
#include <string>
#include <vector>

constexpr int maxpred = 17;
QCMD_EXTERN int nualt;
QCMD_EXTERN int maxualt;
QCMD_EXTERN double gear[maxpred];
QCMD_EXTERN double aspc[maxpred];
QCMD_EXTERN double bpred[maxpred];
QCMD_EXTERN double bpredp[maxpred];
QCMD_EXTERN double bpreds[maxpred];
QCMD_EXTERN double bpredps[maxpred];
QCMD_EXTERN std::vector<std::vector<std::vector<double>>> udalt;
QCMD_EXTERN std::vector<std::vector<std::vector<double>>> upalt;
QCMD_EXTERN std::vector<std::vector<std::vector<double>>> usalt;
QCMD_EXTERN std::vector<std::vector<std::vector<double>>> upsalt;
QCMD_EXTERN bool use_pred;
QCMD_EXTERN std::string polpred;
