// Author: Moses KJ Chung
// Year:   2024

#include "basefile.h"
#include "fatal.h"
#include "files.h"
#include "getcartqm.h"
#include "inform.h"
#include "inquire.h"
#include "nextarg.h"
#include "output.h"
#include "readxyzqm.h"
#include "suffix.h"
#include <iostream>
#include <fstream>
#include <string>

namespace polmdqc
{
///////////////////////////////////////////////////////
//                                                   //
//  getcartqm  --  get a Cartesian coordinates file  //
//                                                   //
///////////////////////////////////////////////////////

// "getcartqm" asks for a Cartesian coordinate file name,
// then reads the formatted or binary coordinates file

void getcartqm(std::ifstream& ffile)
{
    std::string xyzfile;
    bool exist;

    // try to get a filename from the command line arguments
    nextarg(xyzfile, exist);
    if (exist) {
        basefile(xyzfile);
        suffix(xyzfile, "qxyz", "old");
        exist = inquireFile(xyzfile);
    }

    // ask for the user specified input structure filename
    int nask = 0;
    while (!exist and nask<maxask) {
        nask++;
        printf("\n Enter Cartesian Coordinate File Name :  ");
        std::getline(std::cin, xyzfile);
        basefile(xyzfile);
        suffix(xyzfile, "qxyz", "old");
        exist = inquireFile(xyzfile);
    }
    if (!exist) fatal();

    // get file format type by inspection of first character
    // only support archive for getcartqm
    filename = xyzfile;
    coordtype = "CARTESIAN";
    archive = true;
    binary = (!archive);

    // read initial Cartesian coordinates from formatted file
    ffile.open(xyzfile);
    ffile.seekg(0, std::ios::beg);
    readxyzqm(ffile);

    // quit if the Cartesian coordinates file contains no atoms
    if (informAbort) {
        printf("\n GETCARTQM  --  Cartesian Coordinate File was not Read Correctly\n");
        ffile.close();
        fatal();
    }
}
}
