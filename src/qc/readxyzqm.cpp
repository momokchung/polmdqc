// Author: Moses KJ Chung
// Year:   2024

#include "atoms.h"
#include "chkxyz.h"
#include "fatal.h"
#include "files.h"
#include "groupqm.h"
#include "inform.h"
#include "inquire.h"
#include "katoms.h"
#include "readxyzqm.h"
#include "upcase.h"
#include "version.h"
#include <sstream>
#include <vector>

namespace polmdqc
{
//////////////////////////////////////////////////////
//                                                  //
//  readxyzqm  --  input of XYZ-format coordinates  //
//                                                  //
//////////////////////////////////////////////////////

// "readxyzqm" gets a set of Cartesian coordinates from
// an external disk file

void readxyzqm(std::ifstream& ffile)
{
    int size;
    int chg,mult;
    realq xi,yi,zi;
    bool exist,opened,quit;
    std::string xyzfile;
    std::string record;
    std::string sym;
    std::istringstream iss;
    std::vector<realq> xvec,yvec,zvec;
    std::vector<std::string> symvec;
    std::vector<int> grpqvec;
    std::vector<int> chgvec,multvec;

    // initialize the total number of atoms in the system
    n = 0;

    // initialize the total number of groups in the system
    ngrpq = 0;

    // open the input file if it has not already been done
    opened = inquireUnit(ffile);
    if (!opened) {
        xyzfile = filename.substr(0,leng) + ".qxyz";
        version(xyzfile,"old");
        exist = inquireFile(xyzfile);
        if (exist) {
            ffile.open(xyzfile);
        }
        else {
            printf("\n READXYZQM  --  Unable to Find the Cartesian Coordinates File\n");
            fatal();
        }
    }

    // read first line and return if already at end of file
    size = 0;
    while (size == 0) {
        if (std::getline(ffile, record)) {
            iss.clear();
            iss.str(record);
            std::string firstWord = "";
            iss >> firstWord;
            size = firstWord.length();
        }
    }
    
    // parse the first nonempty line to get charge and multiplicity
    iss.clear();
    iss.str(record);
    if (!(iss >> chg >> mult)) {
        printf("\n READXYZQM  --  Error in Cartesian Coordinates File Format\n");
        ffile.close();
        fatal();
    }

    // read the file line by line
    while (std::getline(ffile, record)) {

        // skip empty line
        if (record.empty()) continue;

        // record group information at "--"
        if (record == "--") {
            chgvec.push_back(chg);
            multvec.push_back(mult);
            ngrpq++;
            continue;
        }

        // read in coordinates
        iss.clear();
        iss.str(record);
        if (iss >> sym >> xi >> yi >> zi) {
            upcase(sym);
            symvec.push_back(sym);
            xvec.push_back(xi);
            yvec.push_back(yi);
            zvec.push_back(zi);
            grpqvec.push_back(ngrpq);
            n++;
            continue;
        }

        // read in charge and multiplicity
        iss.clear();
        iss.str(record);
        iss >> chg >> mult;
    }

    // record group information
    chgvec.push_back(chg);
    multvec.push_back(mult);
    ngrpq++;

    // quit if vector sizes don't match
    quit = false;
    if (n != symvec.size()) quit = true;
    if (n != xvec.size()) quit = true;
    if (n != yvec.size()) quit = true;
    if (n != zvec.size()) quit = true;
    if (n != grpqvec.size()) quit = true;
    if (ngrpq != chgvec.size()) quit = true;
    if (ngrpq != multvec.size()) quit = true;
    if (quit) {
        printf("\n READXYZQM  --  Error in Cartesian Coordinates File Format\n");
        ffile.close();
        fatal();
    }

    // allocate global arrays from module atoms
    x.allocate(n);
    y.allocate(n);
    z.allocate(n);

    // allocate global arrays from module katoms
    symbol.allocate(n);

    // allocate global arrays from module groupqm
    grpqlist.allocate(n);
    grpqchg.allocate(ngrpq);
    grpqmult.allocate(ngrpq);

    // assign charge and multiplicity for each group
    for (int i = 0; i < ngrpq; i++) {
        grpqchg[i] = chgvec[i];
        grpqmult[i] = multvec[i];
    }

    // assign symbol and coordinates for each atom
    for (int i = 0; i < n; i++) {
        symbol[i] = symvec[i];
        x[i] = xvec[i];
        y[i] = yvec[i];
        z[i] = zvec[i];
        grpqlist[i] = grpqvec[i];
    }

    // check for atom pairs with identical coordinates
    bool clash = false;
    if (n < 10000) chkxyz(clash);
}
}
