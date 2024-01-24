// Author: Moses KJ Chung
// Year:   2023

#include "active.h"
#include "alphamol.h"
#include "alphmol.h"
#include "atomid.h"
#include "atoms.h"
#include "field.h"
#include "files.h"
#include "final.h"
#include "getcart.h"
#include "inform.h"
#include "initial.h"
#include "katom.h"
#include "kvdw.h"
#include "mechanic.h"
#include "nextarg.h"
#include "readcart.h"
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

// "spacefill" computes the surface area, volume, mean curvature
// gaussian curvature and their respective gradients

void resizeNumDer();

void spacefill(int argc, char** argv)
{
    int mode;
    int frame,nold;
    int numSpaces;
    int width;
    int precision;
    real value,exclude;
    real eps,eps0;
    std::vector<real> denorm;
    std::vector<real> ndenorm;
    std::vector<real> told;
    bool exist,query,first;
    bool doanalyt,donumer;
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
    if (string != "Y") {
        for (int i = 0; i < n; i++) {
            if (atomic[i] == 1) use[i+1] = false;
        }
    }

    // decide whether to compute analytical derivatives
    doanalyt = false;
    nextarg(string,exist);
    if (!exist) {
        printf("\n Compute the Analytical Gradient Vector [N] :  ");
        std::getline(std::cin, string);
    } 
    upcase(string);
    if (string != "Y") {
        doanalyt = true;
    }

    // decide whether to compute numerical derivatives
    donumer = false;
    nextarg(string,exist);
    if (!exist) {
        printf("\n Compute the Numerical Gradient Vector [N] :  ");
        std::getline(std::cin, string);
    } 
    upcase(string);
    if (string != "Y") {
        donumer = true;
    }

    // get the stepsize for numerical gradient calculation
    if (donumer) {
        eps = -1.;
        eps0 = 0.00001;
        query = true;
        nextarg(string,exist);
        if (exist) {
            iss.clear();
            iss.str(string);
            iss >> eps;
            query = false;
        }
        if (query) {
            printf("\n Enter Finite Difference Stepsize [%8.1e Ang] :  ", eps0);
            std::getline(std::cin, record);
            iss.clear();
            iss.str(record);
            iss >> eps;
        }
        if (eps < 0.) eps = eps0;
    }

    // perform dynamic allocation of some local arrays
    denorm.resize(n);
    if (donumer) {
        ndenorm.resize(n);
        resizeNumDer();
    }

    // perform analysis for each successive coordinate structure
    frame = 0;
    while (!informAbort) {
        frame++;
        if (frame > 1) {
            printf("\n Analysis for Archive Structure :        %8d\n", frame);
            if (nold != n) {
                mechanic();
                denorm.resize(n);
                if (donumer) {
                    ndenorm.resize(n);
                    resizeNumDer();
                }
            }
            else {
                for (int i = 0; i < n; i++) {
                    if (type[i] != told[i]) {
                        mechanic();
                        break;                        
                    }
                }
            }
        }

        // compute volume, surface area, mean, and gaussian curvature
        alphamol(exclude, doanalyt);

        if (mode == 1) {
            printf("\n Van der Waals Surface Area and Volume :\n");
        }
        else if (mode == 2) {
            printf("\n Accessible Surface Area and Excluded Volume :\n");
        }

        if (digits >= 8) {
            width = 24;
            precision = 8;
        }
        else if (digits >= 6) {
            width = 22;
            precision = 6;
        }
        else {
            width = 20;
            precision = 4;
        }
        printf("\n Total Area :%*.*f Square Angstroms\n", width, precision, tsurf);

        // attempt to read next structure from the coordinate file
        if (told.size() != n) {
            told.resize(n);
        }
        nold = n;
        for (int i = 0; i < nold; i++) {
            told[i] = type[i];
        }
        first = false;
        readcart(ffile,first);
    }

    // perform any final tasks before program exit
    final();
}

void resizeNumDer() {
    ndsurf.allocate(3*n);
    ndvol.allocate(3*n);
    ndmean.allocate(3*n);
    ndgauss.allocate(3*n);
}
}
