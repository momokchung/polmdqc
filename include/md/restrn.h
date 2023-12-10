// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"
#include <string>
#include <vector>

namespace polmdqc
{
/////////////////////////////////////////////////////////
//                                                     //
//  restrn  --  parameters for geometrical restraints  //
//                                                     //
/////////////////////////////////////////////////////////

// maxfix     maximum number of geometric restraint entries
// npfix      number of position restraints to be applied
// ndfix      number of distance restraints to be applied
// nafix      number of angle restraints to be applied
// ntfix      number of torsional restraints to be applied
// ngfix      number of group distance restraints to be applied
// nchir      number of chirality restraints to be applied
// ipfix      atom number involved in each position restraint
// kpfix      flags to use x-, y-, z-coordinate position restraints
// idfix      atom numbers defining each distance restraint
// iafix      atom numbers defining each angle restraint
// itfix      atom numbers defining each torsional restraint
// igfix      group numbers defining each group distance restraint
// ichir      atom numbers defining each chirality restraint
// depth      depth of shallow Gaussian basin restraint
// width      exponential width coefficient of Gaussian basin
// rflat      flat bottom radius for Gaussian basin restraint
// rwall      radius of spherical droplet boundary restraint
// xpfix      x-coordinate target for each restrained position
// ypfix      y-coordinate target for each restrained position
// zpfix      z-coordinate target for each restrained position
// pfix       force constant and flat-well range for each position
// dfix       force constant and target range for each distance
// afix       force constant and target range for each angle
// tfix       force constant and target range for each torsion
// gfix       force constant and target range for each group distance
// chir       force constant and target range for chiral centers
// use_basin  logical flag governing use of Gaussian basin
// use_wall   logical flag governing use of droplet boundary

MDQC_EXTERN int maxfix;
MDQC_EXTERN int npfix,ndfix;
MDQC_EXTERN int nafix,ntfix;
MDQC_EXTERN int ngfix,nchir;
MDQC_EXTERN std::vector<int> ipfix;
MDQC_EXTERN std::vector<std::vector<int>> kpfix;
MDQC_EXTERN std::vector<std::vector<int>> idfix;
MDQC_EXTERN std::vector<std::vector<int>> iafix;
MDQC_EXTERN std::vector<std::vector<int>> itfix;
MDQC_EXTERN std::vector<std::vector<int>> igfix;
MDQC_EXTERN std::vector<std::vector<int>> ichir;
MDQC_EXTERN double depth,width;
MDQC_EXTERN double rflat,rwall;
MDQC_EXTERN std::vector<double> xpfix;
MDQC_EXTERN std::vector<double> ypfix;
MDQC_EXTERN std::vector<double> zpfix;
MDQC_EXTERN std::vector<std::vector<double>> pfix;
MDQC_EXTERN std::vector<std::vector<double>> dfix;
MDQC_EXTERN std::vector<std::vector<double>> afix;
MDQC_EXTERN std::vector<std::vector<double>> tfix;
MDQC_EXTERN std::vector<std::vector<double>> gfix;
MDQC_EXTERN std::vector<std::vector<double>> chir;
MDQC_EXTERN bool use_basin,use_wall;
}
