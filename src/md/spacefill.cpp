// Author: Moses KJ Chung
// Year:   2023

#include "active.h"
#include "alphamol.h"
#include "atomid.h"
#include "atoms.h"
#include "field.h"
#include "files.h"
#include "getcart.h"
#include "initial.h"
#include "katom.h"
#include "kvdw.h"
#include "nextarg.h"
#include "upcase.h"
#include "usage.h"
#include <iostream>
#include <sstream>

namespace polmdqc
{
///////////////////////////////////////////////////////
//                                                   //
//  spacefill  --  surface area and volume of model  //
//                                                   //
///////////////////////////////////////////////////////

void spacefill(int argc, char** argv)
{
    int mode;
    real value,exclude;
    bool exist,query;
    std::string xyzfile,record,string;
    std::istringstream iss;

    // get the Cartesian coordinates for the system
    initial(argc, argv);
    getcart(ffile);

    // determine the atoms to be used in computation;
    // radii can be changed via the keyword mechanism
    field();
    active();
    katom();
    kvdw();

    query = true;
    nextarg(string,exist);
    if (exist) {
        iss.clear();
        iss.str(string);
        if (iss >> mode) query = false;
    }
    if (query) {
        printf("\n Three Types of Area and Volume can be Computed :");
        printf("\n\n    (1) Van der Waals Area and Volume");
        printf("\n    (2) Accessible Area and Excluded Volume");
        printf("\n\n Enter the Number of your Choice [1] :  ");
        std::getline(std::cin, string);
        iss.clear();
        iss.str(string);
        iss >> mode;
    }
    if (mode!=2) mode = 1;

    // set the excluded/accessible probes
    value = 0.;
    exclude = 0.;
    if (mode == 2) {
        query = true;
        nextarg(string,exist);
        if (exist) {
            iss.clear();
            iss.str(string);
            if (iss >> value) query = false;
        }
        if (query) {
            printf("\n Enter a Value for the Probe Radius [1.4 Ang] :  ");
            std::getline(std::cin, string);
            iss.clear();
            iss.str(string);
            iss >> value;
        }
        if (value < 0.) value = -value;
        if (value == 0.) value = 1.4;
        if (mode == 2) exclude = value;
    }

    // decide whether to include hydrogens in the calculation
    nextarg(string,exist);
    if (!exist) {
        printf("\n Include the Hydrogen Atoms in Computation [N] :  ");
        std::getline(std::cin, string);
    }
    upcase(string);
    if (string!= "Y") {
        for (int i = 0; i < n; i++) {
            if (atomic[i] == 1) use[i+1] = false;
        }
    }

    // compute volume, surface area, mean, and gaussian curvature
    // using AlphaMol
    int flag_deriv = 0;
    alphamol(exclude, flag_deriv);
}
}
