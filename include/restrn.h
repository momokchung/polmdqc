///////////////////////////////////////////////////////////
//                                                       //
//  restrn.h  --  parameters for geometrical restraints  //
//                                                       //
///////////////////////////////////////////////////////////

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


#pragma once
#include "macro.h"
#include <string>
#include <vector>

QCMD_EXTERN int maxfix;
QCMD_EXTERN int npfix,ndfix;
QCMD_EXTERN int nafix,ntfix;
QCMD_EXTERN int ngfix,nchir;
QCMD_EXTERN std::vector<int> ipfix;
QCMD_EXTERN std::vector<std::vector<int>> kpfix;
QCMD_EXTERN std::vector<std::vector<int>> idfix;
QCMD_EXTERN std::vector<std::vector<int>> iafix;
QCMD_EXTERN std::vector<std::vector<int>> itfix;
QCMD_EXTERN std::vector<std::vector<int>> igfix;
QCMD_EXTERN std::vector<std::vector<int>> ichir;
QCMD_EXTERN double depth,width;
QCMD_EXTERN double rflat,rwall;
QCMD_EXTERN std::vector<double> xpfix;
QCMD_EXTERN std::vector<double> ypfix;
QCMD_EXTERN std::vector<double> zpfix;
QCMD_EXTERN std::vector<std::vector<double>> pfix;
QCMD_EXTERN std::vector<std::vector<double>> dfix;
QCMD_EXTERN std::vector<std::vector<double>> afix;
QCMD_EXTERN std::vector<std::vector<double>> tfix;
QCMD_EXTERN std::vector<std::vector<double>> gfix;
QCMD_EXTERN std::vector<std::vector<double>> chir;
QCMD_EXTERN bool use_basin,use_wall;
