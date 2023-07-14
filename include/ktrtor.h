///////////////////////////////////////////////////////////
//                                                       //
//  ktrtor.h  --  torsion-torsion forcefield parameters  //
//                                                       //
///////////////////////////////////////////////////////////

// maxntt    maximum number of torsion-torsion parameter entries
// maxtgrd   maximum dimension of torsion-torsion spline grid
// maxtgrd2  maximum number of torsion-torsion spline grid points
// tnx       number of columns in torsion-torsion spline grid
// tny       number of rows in torsion-torsion spline grid
// ttx       angle values for first torsion of spline grid
// tty       angle values for second torsion of spline grid
// tbf       function values at points on spline grid
// tbx       gradient over first torsion of spline grid
// tby       gradient over second torsion of spline grid
// tbxy      Hessian cross components over spline grid
// ktt       string of torsion-torsion atom classes


#pragma once
#include "macro.h"
#include <string>
#include <vector>

QCMD_EXTERN int maxntt;
constexpr int maxtgrd = 30;
constexpr int maxtgrd2= maxtgrd*maxtgrd;
QCMD_EXTERN std::vector<int> tnx;
QCMD_EXTERN std::vector<int> tny;
QCMD_EXTERN std::vector<std::vector<double>> ttx;
QCMD_EXTERN std::vector<std::vector<double>> tty;
QCMD_EXTERN std::vector<std::vector<double>> tbf;
QCMD_EXTERN std::vector<std::vector<double>> tbx;
QCMD_EXTERN std::vector<std::vector<double>> tby;
QCMD_EXTERN std::vector<std::vector<double>> tbxy;
QCMD_EXTERN std::vector<std::string> ktt;
