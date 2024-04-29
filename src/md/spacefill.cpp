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
#include "initalf.h"
#include "inform.h"
#include "initial.h"
#include "libfunc.h"
#include "katom.h"
#include "kvdw.h"
#include "nextarg.h"
#include "readcart.h"
#include "upcase.h"
#include "usage.h"
#include <iostream>
#include <sstream>
#include <vector>

namespace polmdqc
{
///////////////////////////////////////////////////////
//                                                   //
//  spacefill  --  surface area and volume of model  //
//                                                   //
///////////////////////////////////////////////////////

// "spacefill" computes the surface area, volume, mean curvature
// gaussian curvature and their respective gradients

static void resizeNumDer();

static void printAtomic();

static void printGrad(bool doanalyt, bool donumer, std::vector<real>& denorm, std::vector<real>& ndenorm,
    MDQCArray<real>& der, MDQCArray<real>& nder, std::string& gradstr1, std::string& gradstr2);

void spacefill(int argc, char** argv)
{
    int mode;
    int frame,nold;
    int numSpaces;
    int width;
    int precision;
    real exclude;
    real eps,eps0;
    real old;
    real tsurf0,tvol0,tmean0,tgauss0;
    std::vector<real> denorm;
    std::vector<real> ndenorm;
    std::vector<real> told;
    bool exist,query,first;
    bool dofull,doanalyt,donumer;
    std::string xyzfile,record,string;
    std::istringstream iss;

    // get the Cartesian coordinates for the system
    initial(argc, argv);
    getcart(ffile);

    // determine the atoms to be used in computation;
    // radii can be changed via the keyword mechanism
    active();
    field();
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
        printf("\n Two Types of Area and Volume can be Computed :\n");
        printf("\n    (1) Van der Waals Area and Volume");
        printf("\n    (2) Accessible Area and Excluded Volume");
        printf("\n\n Enter the Number of your Choice [1] :  ");
        std::getline(std::cin, string);
        iss.clear();
        iss.str(string);
        iss >> mode;
    }
    if (mode!=2) mode = 1;

    // set the excluded/accessible probes
    exclude = 0.;
    if (mode == 2) {
        query = true;
        nextarg(string,exist);
        if (exist) {
            iss.clear();
            iss.str(string);
            if (iss >> exclude) query = false;
        }
        if (query) {
            printf("\n Enter a Value for the Probe Radius [1.4 Ang] :  ");
            std::getline(std::cin, string);
            iss.clear();
            iss.str(string);
            iss >> exclude;
        }
        if (exclude < 0.) exclude = -exclude;
        if (exclude == 0.) exclude = 1.4;
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

    // decide whether to output area/volume/curvature for each atom
    dofull = true;
    nextarg(string,exist);
    if (!exist) {
        printf("\n Output Breakdown for Each Atom [N] :  ");
        std::getline(std::cin, string);
    }
    upcase(string);
    if (string != "Y") {
        dofull = false;
    }

    // decide whether to compute analytical derivatives
    doanalyt = true;
    nextarg(string,exist);
    if (!exist) {
        printf("\n Compute the Analytical Gradient Vector [N] :  ");
        std::getline(std::cin, string);
    }
    upcase(string);
    if (string != "Y") {
        doanalyt = false;
    }

    // decide whether to compute numerical derivatives
    donumer = true;
    nextarg(string,exist);
    if (!exist) {
        printf("\n Compute the Numerical Gradient Vector [N] :  ");
        std::getline(std::cin, string);
    }
    upcase(string);
    if (string != "Y") {
        donumer = false;
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

    // initialize AlphaMol
    initalf(1, 1, exclude, doanalyt);

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
            if (!test) printf("\n Analysis for Archive Structure :        %8d\n", frame);
            if (nold != n) {
                active();
                field();
                katom();
                kvdw();
                initalf(1, 1, exclude, doanalyt);
                denorm.resize(n);
                if (donumer) {
                    ndenorm.resize(n);
                    resizeNumDer();
                }
            }
            else {
                for (int i = 0; i < n; i++) {
                    if (type[i] != told[i]) {
                        active();
                        field();
                        katom();
                        kvdw();
                        initalf(1, 1, exclude, doanalyt);
                        break;                        
                    }
                }
            }
        }

        // compute surface area, volume, mean, and gaussian curvature
        alphamol(doanalyt);

        // print atomic surface area, volume, mean, and gaussian curvature
        if (dofull and !test) {
            printAtomic();
        }

        // print total surface area, volume, mean, and gaussian curvature
        if (!test) {
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

            printf("\n Total Area :              %*.*f Square Angstroms", width, precision, tsurf);
            printf("\n Total Volume :            %*.*f Square Angstroms", width, precision, tvol);
            printf("\n Total Mean Curvature :    %*.*f Square Angstroms", width, precision, tmean);
            printf("\n Total Gaussian Curvature :%*.*f Square Angstroms\n", width, precision, tgauss);
        }

        // get the Cartesian component two-sided numerical gradients
        for (int i = 0; i < n; i++) {
            if (donumer and use[i+1]) {
                for (int j = 0; j < 3; j++) {
                    if (j == 0) {
                        old = x[i];
                        x[i] -= (real)0.5 * eps;
                    }
                    else if (j == 1) {
                        old = y[i];
                        y[i] -= (real)0.5 * eps;
                    }
                    else if (j == 2) {
                        old = z[i];
                        z[i] -= (real)0.5 * eps;
                    }
                    alphamol(false);
                    tsurf0 = tsurf;
                    tvol0 = tvol;
                    tmean0 = tmean;
                    tgauss0 = tgauss;
                    if (j == 0) {
                        x[i] += eps;
                    }
                    else if (j == 1) {
                        y[i] += eps;
                    }
                    else if (j == 2) {
                        z[i] += eps;
                    }
                    alphamol(false);
                    if (j == 0) {
                        x[i] = old;
                    }
                    else if (j == 1) {
                        y[i] = old;
                    }
                    else if (j == 2) {
                        z[i] = old;
                    }
                    ndsurf[3*i+j] = (tsurf - tsurf0) / eps;
                    ndvol[3*i+j] = (tvol - tvol0) / eps;
                    ndmean[3*i+j] = (tmean - tmean0) / eps;
                    ndgauss[3*i+j] = (tgauss - tgauss0) / eps;
                }
            }
        }

        // print the total gradient components for each atom
        if ((doanalyt or donumer) and !test) {
            std::string gradstr1,gradstr2;
            gradstr1 = "Surface Area";
            gradstr2 = "SA";
            printGrad(doanalyt, donumer, denorm, ndenorm, dsurf, ndsurf, gradstr1, gradstr2);

            gradstr1 = "Volume";
            gradstr2 = "Vol";
            printGrad(doanalyt, donumer, denorm, ndenorm, dvol, ndvol, gradstr1, gradstr2);

            gradstr1 = "Mean Curvature";
            gradstr2 = "MC";
            printGrad(doanalyt, donumer, denorm, ndenorm, dmean, ndmean, gradstr1, gradstr2);

            gradstr1 = "Gaussian Curvature";
            gradstr2 = "GC";
            printGrad(doanalyt, donumer, denorm, ndenorm, dgauss, ndgauss, gradstr1, gradstr2);
        }

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
    ffile.close();
    if (!test) final();
}

static void resizeNumDer() {
    ndsurf.allocate(3*n);
    ndvol.allocate(3*n);
    ndmean.allocate(3*n);
    ndgauss.allocate(3*n);
}

static void printAtomic()
{
    printf("\n Surface Area, Volume, and Mean/Guassian Curvature of Individual Atoms :\n");
    int s1;
    if (digits >= 8) {
        s1 = 12;
    }
    else if (digits >= 6) {
        s1 = 10;
    }
    else {
        s1 = 8;
    }
    printf("\n     Atom%*sS Area%*sVolume%*sM Curv%*sG Curv\n\n", s1,"",s1,"",s1,"",s1,"");
    for (int i = 0; i < n; i++) {
        if (digits >= 8) {
                printf(" %8d%18.8f%18.8f%18.8f%18.8f\n", i+1,surf[i],vol[i],mean[i],gauss[i]);
            }
            else if (digits >= 6) {
                printf(" %8d%16.6f%16.6f%16.6f%16.6f\n", i+1,surf[i],vol[i],mean[i],gauss[i]);
            }
            else {
                printf(" %8d%14.4f%14.4f%14.4f%14.4f\n", i+1,surf[i],vol[i],mean[i],gauss[i]);
            }
    }
}

static void printGrad(bool doanalyt, bool donumer, std::vector<real>& denorm, std::vector<real>& ndenorm,
    MDQCArray<real>& der, MDQCArray<real>& nder, std::string& gradstr1, std::string& gradstr2)
{
    real totnorm,ntotnorm,rms,nrms;

    if (digits >= 8) {
        printf("\n----------------------------------------------------------------------------");
    }
    else if (digits >= 6) {
        printf("\n---------------------------------------------------------------------------");
    }
    else {
        printf("\n-------------------------------------------------------------------------");
    }
    printf("\n %s Gradient over Individual Atoms :\n", gradstr1.c_str());
    int s1,s2,s3,s4,s5;
    if (digits >= 8) {
        s1 = 4; s2 = 12; s3 = 11; s4 = 11; s5 = 12;
    }
    else if (digits >= 6) {
        s1 = 6; s2 = 12; s3 = 9; s4 = 9; s5 = 12;
    }
    else {
        s1 = 6; s2 = 14; s3 = 7; s4 = 7; s5 = 10;
    }
    printf("\n  Type%*sAtom%*sdE/dX%*sdE/dY%*sdE/dZ%*sNorm\n\n", s1,"",s2,"",s3,"",s4,"",s5,"");

    totnorm = 0.;
    ntotnorm = 0.;
    for (int i = 0; i < n; i++) {
        if (doanalyt) {
            denorm[i] = REAL_POW(der[3*i+0],2) + REAL_POW(der[3*i+1],2) + REAL_POW(der[3*i+2],2);
            totnorm = totnorm + denorm[i];
            denorm[i] = REAL_SQRT(denorm[i]);
            if (digits >= 8) {
                printf(" Anlyt%8d %16.8f%16.8f%16.8f%16.8f\n", i+1,der[3*i+0],der[3*i+1],der[3*i+2],denorm[i]);
            }
            else if (digits >= 6) {
                printf(" Anlyt  %8d   %14.6f%14.6f%14.6f  %14.6f\n", i+1,der[3*i+0],der[3*i+1],der[3*i+2],denorm[i]);
            }
            else {
                printf(" Anlyt  %8d       %12.4f%12.4f%12.4f  %12.4f\n", i+1,der[3*i+0],der[3*i+1],der[3*i+2],denorm[i]);
            }
        }
        if (donumer) {
            ndenorm[i] = REAL_POW(nder[3*i+0],2) + REAL_POW(nder[3*i+1],2) + REAL_POW(nder[3*i+2],2);
            ntotnorm = ntotnorm + ndenorm[i];
            ndenorm[i] = REAL_SQRT(ndenorm[i]);
            if (digits >= 8) {
                printf(" Numer%8d %16.8f%16.8f%16.8f%16.8f\n", i+1,nder[3*i+0],nder[3*i+1],nder[3*i+2],ndenorm[i]);
            }
            else if (digits >= 6) {
                printf(" Numer  %8d   %14.6f%14.6f%14.6f  %14.6f\n", i+1,nder[3*i+0],nder[3*i+1],nder[3*i+2],ndenorm[i]);
            }
            else {
                printf(" Numer  %8d       %12.4f%12.4f%12.4f  %12.4f\n", i+1,nder[3*i+0],nder[3*i+1],nder[3*i+2],ndenorm[i]);
            }
        }
    }

    // print the total norm for the analytical gradient
    if (doanalyt or donumer) {
        printf("\n %s Gradient Norm and RMS per Atom :\n\n", gradstr1.c_str());
    }
    if (doanalyt) {
        totnorm = REAL_SQRT(totnorm);
        if (digits >= 8) {
            printf(" Anlyt      %s Gradient Norm Value            %20.8f\n", gradstr2.c_str(), totnorm);
        }
        else if (digits >= 6) {
            printf(" Anlyt      %s Gradient Norm Value            %18.6f\n", gradstr2.c_str(), totnorm);
        }
        else {
            printf(" Anlyt      %s Gradient Norm Value            %16.4f\n", gradstr2.c_str(), totnorm);
        }
    }

    // print the total norm for the numerical gradient
    if (donumer) {
        ntotnorm = REAL_SQRT(ntotnorm);
        if (digits >= 8) {
            printf(" Numer      %s Gradient Norm Value            %20.8f\n", gradstr2.c_str(), ntotnorm);
        }
        else if (digits >= 6) {
            printf(" Numer      %s Gradient Norm Value            %18.6f\n", gradstr2.c_str(), ntotnorm);
        }
        else {
            printf(" Numer      %s Gradient Norm Value            %16.4f\n", gradstr2.c_str(), ntotnorm);
        }
    }

    printf("\n");

    // print the rms per atom norm for the analytical gradient
    if (doanalyt) {
        rms = totnorm / REAL_SQRT(static_cast<real>(nuse));
        if (digits >= 8) {
            printf(" Anlyt      RMS %s Gradient over All Atoms    %20.8f\n", gradstr2.c_str(), rms);
        }
        else if (digits >= 6) {
            printf(" Anlyt      RMS %s Gradient over All Atoms    %18.6f\n", gradstr2.c_str(), rms);
        }
        else {
            printf(" Anlyt      RMS %s Gradient over All Atoms    %16.4f\n", gradstr2.c_str(), rms);
        }
    }

    // print the rms per atom norm for the numerical gradient
    if (donumer) {
        nrms = ntotnorm / REAL_SQRT(static_cast<real>(nuse));
        if (digits >= 8) {
            printf(" Numer      RMS %s Gradient over All Atoms    %20.8f\n", gradstr2.c_str(), nrms);
        }
        else if (digits >= 6) {
            printf(" Numer      RMS %s Gradient over All Atoms    %18.6f\n", gradstr2.c_str(), nrms);
        }
        else {
            printf(" Numer      RMS %s Gradient over All Atoms    %16.4f\n", gradstr2.c_str(), nrms);
        }
    }
}
}
