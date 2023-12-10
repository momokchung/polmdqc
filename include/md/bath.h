// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"
#include <string>

namespace polmdqc
{
////////////////////////////////////////////////////////
//                                                    //
//  bath  --  thermostat and barostat control values  //
//                                                    //
////////////////////////////////////////////////////////

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

constexpr int maxnose = 4;
MDQC_EXTERN int voltrial;
MDQC_EXTERN double kelvin,atmsph;
MDQC_EXTERN double tautemp,taupres;
MDQC_EXTERN double compress,collide;
MDQC_EXTERN double eta,volmove;
MDQC_EXTERN double vbar,qbar,gbar;
MDQC_EXTERN double vnh[maxnose];
MDQC_EXTERN double qnh[maxnose];
MDQC_EXTERN double gnh[maxnose];
MDQC_EXTERN bool isothermal;
MDQC_EXTERN bool isobaric;
MDQC_EXTERN bool anisotrop;
MDQC_EXTERN std::string volscale;
MDQC_EXTERN std::string barostat;
MDQC_EXTERN std::string thermostat;
}
