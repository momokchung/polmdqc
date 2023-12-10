// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include "macro.h"

namespace polmdqc
{
//////////////////////////////////////////////////////////
//                                                      //
//  units  --  physical constants and unit conversions  //
//                                                      //
//////////////////////////////////////////////////////////

// literature references:

// M. Stock, R. Davis, E. de Mirandes and M. J. T. Milton, "The
// Revision of the SI - The Result of Three Decades of Progress
// in Metrology", Metrologia, 56, 022001 (2019)

// E. Teisinga, P. J. Mohr, D. B. Newell and B. N. Taylor, "CODATA
// Recommended Values of the Fundamental Physical Constants: 2018",
// Journal of Physical and Chemical Reference Data, 50, 033105 (2021)

// Where appropriate, values are from the November 2018 revision
// of SI units to fixed values by the 26th General Conference on
// Weights and Measures; other values are taken from 2018 CODATA
// reference constants or are described below

// The conversion from calorie to Joule is the definition of the
// thermochemical calorie as 1 cal = 4.1840 J from ISO 31-4 (1992)

// The "coulomb" energy conversion factor is found by dimensional
// analysis of Coulomb's Law, that is by dividing the square of the
// elementary charge in Coulombs by 4*pi*eps0*rij, where eps0 is
// the permittivity of vacuum (the "electric constant"); note that
// eps0 is typically given in F/m, equivalent to C**2/(J-m)

// The approximate value used for the Debye, 3.33564 x 10-30 C-m,
// is from IUPAC Compendium of Chemical Technology, 2nd Ed. (1997)

// The value of "prescon" is based on definition of 1 atmosphere
// as 101325 Pa set by the 10th Conference Generale des Poids et
// Mesures (Paris, 1954), where a Pascal (Pa) is equal to a J/m**3

// avogadro    Avogadro's number (N) in particles/mole
// lightspd    speed of light in vacuum (c) in cm/ps
// boltzmann   Boltzmann constant (kB) in g*Ang**2/ps**2/mole/K
// gasconst    ideal gas constant (R) in kcal/mole/K
// elemchg     elementary charge of a proton in Coulombs
// vacperm     vacuum permittivity (electric constant, eps0) in F/m
// emass       mass of an electron in atomic mass units
// planck      Planck's constant (h) in J-s
// joule       conversion from calorie to joule
// ekcal       conversion from kcal to g*Ang**2/ps**2
// bohr        conversion from Bohr to Angstrom
// hartree     conversion from Hartree to kcal/mole
// evolt       conversion from Hartree to electron-volt
// efreq       conversion from Hartree to cm-1
// coulomb     conversion from electron**2/Ang to kcal/mole
// elefield    conversion from electron**2/Ang to megavolt/cm
// debye       conversion from electron-Ang to Debye
// prescon     conversion from kcal/mole/Ang**3 to Atm

constexpr double avogadro = 6.02214076e+23;
constexpr double lightspd = 2.99792458e-2;
constexpr double boltzmann = 0.8314462618e0;
constexpr double gasconst = 1.9872042586e-3;
constexpr double elemchg = 1.602176634e-19;
constexpr double vacperm = 8.8541878128e-12;
constexpr double emass = 5.48579909065e-4;
constexpr double planck = 6.62607015e-34;
constexpr double joule = 4.1840e0;
constexpr double ekcal = 4.1840e+2;
constexpr double bohr = 0.529177210903e0;
constexpr double hartree = 627.509474063e0;
constexpr double evolt = 27.211386245988e0;
constexpr double efreq = 2.194746314e+5;
constexpr double coulomb = 332.0637133e0;
constexpr double elefield = 1439.96455e0;
constexpr double debye = 4.80321e0;
constexpr double prescon = 6.85684112e+4;
}
