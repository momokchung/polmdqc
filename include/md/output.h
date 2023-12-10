// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"
#include <string>

namespace polmdqc
{
/////////////////////////////////////////////////////////
//                                                     //
//  output  --  output file format control parameters  //
//                                                     //
/////////////////////////////////////////////////////////

// archive     logical flag for coordinates in Tinker XYZ format
// binary      logical flag for coordinates in DCD binary format
// noversion   logical flag governing use of filename versions
// overwrite   logical flag to overwrite intermediate files inplace
// arcsave     logical flag to save coordinates in Tinker XYZ format
// dcdsave     logical flag to save coordinates in DCD binary format
// cyclesave   logical flag to mark use of numbered cycle files
// velsave     logical flag to save velocity vector components
// frcsave     logical flag to save force vector components
// uindsave    logical flag to save induced atomic dipoles
// coordtype   selects Cartesian, internal, rigid body or none

MDQC_EXTERN bool archive;
MDQC_EXTERN bool binary;
MDQC_EXTERN bool noversion;
MDQC_EXTERN bool overwrite;
MDQC_EXTERN bool cyclesave;
MDQC_EXTERN bool arcsave;
MDQC_EXTERN bool dcdsave;
MDQC_EXTERN bool velsave;
MDQC_EXTERN bool frcsave;
MDQC_EXTERN bool uindsave;
MDQC_EXTERN std::string coordtype;
}
