///////////////////////////////////////////////////////////
//                                                       //
//  output.h  --  output file format control parameters  //
//                                                       //
///////////////////////////////////////////////////////////

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


#pragma once
#include "macro.h"
#include <string>

QCMD_EXTERN bool archive;
QCMD_EXTERN bool binary;
QCMD_EXTERN bool noversion;
QCMD_EXTERN bool overwrite;
QCMD_EXTERN bool cyclesave;
QCMD_EXTERN bool arcsave;
QCMD_EXTERN bool dcdsave;
QCMD_EXTERN bool velsave;
QCMD_EXTERN bool frcsave;
QCMD_EXTERN bool uindsave;
QCMD_EXTERN std::string coordtype;
