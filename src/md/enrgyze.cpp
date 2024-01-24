// Author: Moses KJ Chung
// Year:   2023

#include "atoms.h"
#include "calcMode.h"
#include "energi.h"
#include "energy.h"
#include "enrgyze.h"
#include "inform.h"
#include "inter.h"
#include "mdqclimits.h"
#include "molcul.h"

namespace polmdqc
{
/////////////////////////////////////////////////////
//                                                 //
//  enrgyze  --  compute & report energy analysis  //
//                                                 //
/////////////////////////////////////////////////////

// "enrgyze" is an auxiliary routine for the analyze program
// that performs the energy analysis and prints the total and
// intermolecular energies

void enrgyze(bool test)
{
    constexpr CalcMode CalculationMode = CalcMode::Analysis;

    // perform the energy analysis by atom and component
    energy<CalculationMode>();

    // intermolecular energy for systems with multiple molecules
    int numSpaces = 9;
    int width = 16;
    int precision = 4;
    bool useScientific = false;
    if (digits >= 6) {
        numSpaces = 7;
        width = 18;
        precision = 6;
    }
    if (digits >= 8) {
        numSpaces = 5;
        width = 20;
        precision = 8;
    }
    if (std::abs(einter) >= 1.e10) useScientific = true;

    if (nmol>1 and nmol<n and !use_ewald and !test) {
        if (!useScientific) {
            printf("\n Intermolecular Energy :%*s%*.*f Kcal/mole\n", numSpaces, "", width, precision, einter);
        }
        else {
            printf("\n Intermolecular Energy :%*s%*.*e Kcal/mole\n", numSpaces, "", width, precision, einter);
        }
    }

    // print out the total potential energy of the system
    numSpaces = 8;
    width = 16;
    precision = 4;
    useScientific = false;
    if (digits >= 6) {
        numSpaces = 6;
        width = 18;
        precision = 6;
    }
    if (digits >= 8) {
        numSpaces = 4;
        width = 20;
        precision = 8;
    }
    if (std::abs(esum) >= 1.e10) useScientific = true;
    if (!test) {
        if (!useScientific) {
            printf("\n Total Potential Energy :%*s%*.*f Kcal/mole\n", numSpaces, "", width, precision, esum);
        }
        else {
            printf("\n Total Potential Energy :%*s%*.*e Kcal/mole\n", numSpaces, "", width, precision, esum);
        }
    }
}
}
