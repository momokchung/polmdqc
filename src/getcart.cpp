/////////////////////////////////////////////////////////
//                                                     //
//  getcart.cpp  --  get a Cartesian coordinates file  //
//                                                     //
/////////////////////////////////////////////////////////

// "getcart" asks for a Cartesian coordinate file name, then
// reads the formatted or binary coordinates file


#include "basefile.h"
#include "fatal.h"
#include "files.h"
#include "getcart.h"
#include "inform.h"
#include "inquire.h"
#include "nextarg.h"
#include "output.h"
#include "readxyz.h"
#include "suffix.h"
#include <iostream>
#include <fstream>
#include <string>

void getcart(int &ixyz)
{
    std::string xyzfile;
    bool exist;
    // try to get a filename from the command line arguments
    nextarg (xyzfile, exist);
    if (exist) {
        basefile(xyzfile);
        suffix(xyzfile, "xyz", "old");
        exist = inquire(xyzfile);
    }

    // ask for the user specified input structure filename
    int nask = 0;
    while (!exist and nask<maxask) {
        nask++;
        printf("\n Enter Cartesian Coordinate File Name :  ");
        std::getline(std::cin, xyzfile);
        basefile(xyzfile);
        suffix(xyzfile, "xyz", "old");
        exist = inquire(xyzfile);
    }
    if (!exist) fatal();

    // get file format type by inspection of first character
    filename = xyzfile;
    coordtype = "CARTESIAN";
    std::ifstream file(xyzfile);
    char letter = ' ';
    file.get(letter);
    file.close();
    archive = false;
    if (std::isspace(letter)) archive = true;
    if (letter>='0' and letter <='9') archive = true;
    binary = (!archive);

    // read initial Cartesian coordinates from formatted file
    if (archive) {
        readxyz(xyzfile);
    }

    // TODO: read in binary files

    // quit if the Cartesian coordinates file contains no atoms
    if (informAbort) {
        printf("\n GETCART  --  Cartesian Coordinate File was not Read Correctly");
        fatal();
    }
}
