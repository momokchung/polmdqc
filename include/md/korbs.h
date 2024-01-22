// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"
#include <string>

namespace polmdqc
{
/////////////////////////////////////////////////////////
//                                                     //
//  korbs  --  pisystem orbital forcefield parameters  //
//                                                     //
/////////////////////////////////////////////////////////

// maxnpi     maximum number of pisystem bond parameter entries
// maxnpi5    maximum number of 5-membered ring pibond entries
// maxnpi4    maximum number of 4-membered ring pibond entries
// electron   number of pi-electrons for each atom class
// ionize     ionization potential for each atom class
// repulse    repulsion integral value for each atom class
// sslope     slope for bond stretch vs. pi-bond order
// sslope5    slope for 5-ring bond stretch vs. pi-bond order
// sslope4    slope for 4-ring bond stretch vs. pi-bond order
// tslope     slope for 2-fold torsion vs. pi-bond order
// tslope5    slope for 5-ring 2-fold torsion vs. pi-bond order
// tslope4    slope for 4-ring 2-fold torsion vs. pi-bond order
// kpi        string of atom classes for pisystem bonds
// kpi5       string of atom classes for 5-ring pisystem bonds
// kpi4       string of atom classes for 4-ring pisystem bonds

MDQC_EXTERN int maxnpi;
MDQC_EXTERN int maxnpi5;
MDQC_EXTERN int maxnpi4;
MDQC_EXTERN MDQCArray<real> electron;
MDQC_EXTERN MDQCArray<real> ionize;
MDQC_EXTERN MDQCArray<real> repulse;
MDQC_EXTERN MDQCArray<real> sslope;
MDQC_EXTERN MDQCArray<real> sslope5;
MDQC_EXTERN MDQCArray<real> sslope4;
MDQC_EXTERN MDQCArray<real> tslope;
MDQC_EXTERN MDQCArray<real> tslope5;
MDQC_EXTERN MDQCArray<real> tslope4;
MDQC_EXTERN MDQCArray<std::string> kpi;
MDQC_EXTERN MDQCArray<std::string> kpi5;
MDQC_EXTERN MDQCArray<std::string> kpi4;
}
