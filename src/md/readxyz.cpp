// Author: Moses KJ Chung
// Year:   2023

#include "atomid.h"
#include "atoms.h"
#include "boxes.h"
#include "chkxyz.h"
#include "couple.h"
#include "darray.h"
#include "fatal.h"
#include "files.h"
#include "getline.h"
#include "gettext.h"
#include "getword.h"
#include "inform.h"
#include "inquire.h"
#include "lattice.h"
#include "readxyz.h"
#include "titles.h"
#include "trimtext.h"
#include "unitcell.h"
#include "version.h"
#include <algorithm>
#include <sstream>
#include <vector>

namespace polmdqc
{
////////////////////////////////////////////////////
//                                                //
//  readxyz  --  input of XYZ-format coordinates  //
//                                                //
////////////////////////////////////////////////////

// "readxyz" gets a set of Cartesian coordinates from
// an external disk file

void readxyz(std::ifstream& ffile)
{
    int nmax;
    int next,size;
    int failAtom;
    std::vector<int> list;
    bool exist,opened;
    bool quit,reorder;
    bool clash;
    std::string xyzfile;
    std::string record;
    std::string string;
    std::istringstream iss;
    std::istringstream isstmp;

    // initialize the total number of atoms in the system
    n = 0;

    // open the input file if it has not already been done
    opened = inquireUnit(ffile);
    if (!opened) {
        xyzfile = filename.substr(0,leng) + ".xyz";
        version(xyzfile,"old");
        exist = inquireFile(xyzfile);
        if (exist) {
            ffile.open(xyzfile);
        }
        else {
            printf("\n READXYZ  --  Unable to Find the Cartesian Coordinates File\n");
            fatal();
        }
    }

    // read first line and return if already at end of file
    quit = false;
    informAbort = true;
    size = 0;
    while (size == 0) {
        if (std::getline(ffile, record)) {
            iss.clear();
            iss.str(record);
            std::string firstWord = "";
            iss >> firstWord;
            size = firstWord.length();
        }
        else {
            goto label_80;
        }
    }
    informAbort = false;
    quit = true;

    // parse the title line to get the number of atoms
    next = 0;
    gettext(record, string, next);
    iss.clear();
    iss.str(string);
    if (!(iss >> n)) goto label_80;

    // extract the title and determine its length
    string = record.substr(next);
    getline(string);
    ltitle = string.length();
    title = string;
    if (ltitle == 0) title = "";

    // check for too few or too many total atoms in the file
    if (n <= 0) {
        printf("\n READXYZ  --  The Coordinate File Does Not Contain Any Atoms\n");
        fatal();
    }
    else if (n > maxatm) {
        printf("\n READXYZ  --  The Maximum of%9d Atoms has been Exceeded\n", maxatm);
        fatal();
    }

    // allocate global arrays from module atomid
    tag.allocate(n);
    name.allocate(n);

    // allocate global arrays from module atoms
    type.allocate(n);
    x.allocate(n);
    y.allocate(n);
    z.allocate(n);

    // allocate global arrays from module couple
    n12.allocate(n);
    i12.allocate(n);

    // initialize coordinates and connectivities for each atom
    for (int i = 0; i < n; i++) {
        tag[i] = -1;
        name[i] = "   ";
        x[i] = 0.;
        y[i] = 0.;
        z[i] = 0.;
        type[i] = -1;
        n12[i] = 0;
        for (int j = 0; j < maxval; j++) {
            i12[i][j] = -1;
        }
    }

    // read the coordinates and connectivities for each atom
    failAtom = 0;
    for (int i = 0; i < n; i++) {
        failAtom = i+1;
        size = -1;
        while (size == -1) {
            unitcell();
            if (!std::getline(ffile, record)) goto label_80;
            size = trimtext(record);
            if (i == 0) {
                next = 0;
                getword (record, name[i], next);
                if (name[i].length() == 0) {
                    iss.clear();
                    iss.str(record);
                    if ((iss >> xbox >> ybox >> zbox >> alphaA >> betaA >> gammaA)) size = -1;
                }
            }
            lattice();
        }
        iss.clear();
        iss.str(record);
        int tagInt;
        if (iss >> tagInt) tag[i] = tagInt - 1;
        next = 0;
        getword(record, name[i], next);
        string = record.substr(next);
        iss.clear();
        iss.str(string);
        iss >> x[i] >> y[i] >> z[i] >> type[i];
        type[i]--;
        for (int j = 0; j < maxval; j++) {
            int i12Int;
            if (iss >> i12Int) {
                i12[i][j] = i12Int - 1;
            }
            else {
                break;
            }
        }
    }
    quit = false;
    label_80:
    if (!opened) ffile.close();

    // an error occurred in reading the coordinate file
    if (quit) {
        printf("\n READXYZ  --  Error in Coordinate File at Atom%9d\n", failAtom);
        fatal();
    }

    // for each atom, count and sort its attached atoms
    for (int i = 0; i < n; i++) {
        for (int j = maxval-1; j >= 0; j--) {
            if (i12[i][j] != -1) {
                n12[i] = j + 1;
                break;
            }
        }
        std::sort(i12[i], i12[i] + n12[i]);
    }

    // perform dynamic allocation of some local arrays
    nmax = 0;
    for (int i = 0; i < n; i++) {
        nmax = std::max(tag[i]+1, nmax);
        for (int j = 0; j < n12[i]; j++) {
            nmax = std::max(i12[i][j]+1, nmax);
        }
    }
    list.resize(nmax);

    // check for scrambled atom order and attempt to renumber
    reorder = false;
    for (int i = 0; i < n; i++) {
        list[tag[i]] = i;
        if (tag[i] != i) reorder = true;
    }
    if (reorder) {
        if (!test) printf("\n READXYZ  --  Atom Labels not Sequential, Attempting to Renumber\n");
        for (int i = 0; i < n; i++) {
            tag[i] = i;
            for (int j = 0; j < n12[i]; j++) {
                i12[i][j] = list[i12[i][j]];
            }
            std::sort(i12[i], i12[i] + n12[i]);
        }
    }

    // check for atom pairs with identical coordinates
    clash = false;
    if (n < 10000) chkxyz(clash);

    // make sure all atom connectivities are bidirectional
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n12[i]; j++) {
            int k = i12[i][j];
            for (int m = 0; m < n12[k]; m++) {
                if (i12[k][m] == i) goto label_130;
            }
            printf("\n READXYZ  --  Check Connection of Atoms%9d and%9d\n", k+1, i+1);
            fatal();
            label_130:;
        }
    }
}
}
