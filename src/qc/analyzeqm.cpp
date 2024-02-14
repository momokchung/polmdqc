// Author: Moses KJ Chung
// Year:   2024

#include "analyzeqm.h"
#include "files.h"
#include "final.h"
#include "getcartqm.h"
#include "inform.h"
#include "initialqm.h"
#include "mechanicqm.h"
#include "nextarg.h"
#include "trimtext.h"
#include "upcase.h"
#include <iostream>
#include <string>

namespace polmdqc
{
///////////////////////////////////////////////////////
//                                                   //
//  analyzeqm  --  energy partitioning and analysis  //
//                                                   //
///////////////////////////////////////////////////////

// "analyzeqm" computes the quantum mechanical energy

void analyzeqm(int argc, char** argv)
{
    bool exist;
    bool doenergy,dodetail;
    bool domoment,dovirial;
    char letter;
    std::string string;

    // set up the structure and mechanics calculation
    initialqm(argc, argv);
    getcartqm(ffile);
    mechanicqm();

    // get the desired types of analysis to be performed
    nextarg(string,exist);
    if (!exist and !test) {
        printf("\n The PolMDQC QM Energy Analysis Utility Can :\n");
        // printf("\n General System and Force Field Information [G]");
        // printf("\n Force Field Parameters for Interactions [P]");
        printf("\n Total Potential Energy Calculation [E]");
        // printf("\n Energy Breakdown over Each of the Atoms [A]");
        // printf("\n List of the Large Individual Interactions [L]");
        // printf("\n Details for All Individual Interactions [D]");
        // printf("\n Electrostatic Moments and Principle Axes [M]");
        // printf("\n Internal Virial & Instantaneous Pressure [V]");
        // printf("\n Connectivity Lists for Each of the Atoms [C]");
        printf("\n\n Enter the Desired Analysis [E] :  ");
        std::getline(std::cin, string);
    }

    // set option control flags based desired analysis types
    doenergy = false;
    dodetail = false;
    domoment = false;
    dovirial = false;
    upcase(string);
    int stringlen = trimtext(string);
    for (int i = 0; i <= stringlen; i++) {
        letter = string[i];
        if (letter == 'E') doenergy = true;
        if (letter == 'D') dodetail = true;
        if (letter == 'M') domoment = true;
        if (letter == 'V') dovirial = true;
    }

    // setup to write out all of the individual energy terms
    if (dodetail) {
        doenergy = true;
        verbose = true;
        debug = true;
    }
    else {
        debug = false;
    }

    // decide whether to perform analysis of individual frames
    informAbort = true;
    if (doenergy or domoment or dovirial) informAbort = false;

    // // make the call to compute the potential energy
    // if (doenergy) energyqm();

    // perform any final tasks before program exit
    ffile.close();
    if (!test) final();
}
}
