//////////////////////////////////////////////////////////
//                                                      //
//  bath.h  --  thermostat and barostat control values  //
//                                                      //
//////////////////////////////////////////////////////////

// maxnose     maximum length of Nose-Hoover thermostat chain
// voltrial    mean number of steps between Monte Carlo moves
// kelvin      target value for the system temperature (K)
// atmsph      target value for the system pressure (atm)
// tautemp     time constant for Berendsen thermostat (psec)
// taupres     time constant for Berendsen barostat (psec)
// compress    isothermal compressibility of medium (atm-1)
// collide     collision frequency for Andersen thermostat
// eta         velocity value for Bussi-Parrinello barostat
// volmove     maximum volume move for Monte Carlo barostat (Ang**3)
// vbar        velocity of log volume for Nose-Hoover barostat
// qbar        mass of the volume for Nose-Hoover barostat
// gbar        force for the volume for Nose-Hoover barostat
// vnh         velocity of each chained Nose-Hoover thermostat
// qnh         mass for each chained Nose-Hoover thermostat
// gnh         force for each chained Nose-Hoover thermostat
// isothermal  logical flag governing use of temperature control
// isobaric    logical flag governing use of pressure control
// anisotrop   logical flag governing use of anisotropic pressure
// thermostat  choice of temperature control method to be used
// barostat    choice of pressure control method to be used
// volscale    choice of scaling method for Monte Carlo barostat


#pragma once
#include "macro.h"
#include <string>

constexpr int maxnose = 4;
QCMD_EXTERN int voltrial;
QCMD_EXTERN double kelvin,atmsph;
QCMD_EXTERN double tautemp,taupres;
QCMD_EXTERN double compress,collide;
QCMD_EXTERN double eta,volmove;
QCMD_EXTERN double vbar,qbar,gbar;
QCMD_EXTERN double vnh[maxnose];
QCMD_EXTERN double qnh[maxnose];
QCMD_EXTERN double gnh[maxnose];
QCMD_EXTERN bool isothermal;
QCMD_EXTERN bool isobaric;
QCMD_EXTERN bool anisotrop;
QCMD_EXTERN std::string volscale;
QCMD_EXTERN std::string barostat;
QCMD_EXTERN std::string thermostat;
