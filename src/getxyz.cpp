///////////////////////////////////////////////////////
//                                                   //
//  getxyz.cpp  --  get XYZ-format coordinates file  //
//                                                   //
///////////////////////////////////////////////////////

// "getxyz" asks for a Cartesian coordinate file name,
// then reads in the coordinates file


#include "basefile.h"
#include "fatal.h"
#include "files.h"
#include "getxyz.h"
#include "inform.h"
#include "inquire.h"
#include "nextarg.h"
#include "output.h"
#include "readxyz.h"
#include "suffix.h"
#include <string>

void getxyz()
{
    bool exist;
    std::string xyzfile;

    // try to get a filename from the command line arguments
    nextarg (xyzfile, exist);
    if (exist) {
        basefile(xyzfile);
        suffix(xyzfile, "xyz", "old");
        exist = inquireFile(xyzfile);
    }

    // ask for the user specified input structure filename
    int nask = 0;
    while (!exist and nask<maxask) {
        nask++;
        printf("\n Enter Cartesian Coordinate File Name :  ");
        std::getline(std::cin, xyzfile);
        basefile(xyzfile);
        suffix(xyzfile, "xyz", "old");
        exist = inquireFile(xyzfile);
    }
    if (!exist) fatal();

    // first open and then read the Cartesian coordinates file
    filename = xyzfile;
    coordtype = "CARTESIAN";
    std::ifstream ffile(xyzfile);
    readxyz(ffile);
    ffile.close();

    // quit if the Cartesian coordinates file contains no atoms
    if (informAbort) {
        printf("\n GETXYZ  --  Cartesian Coordinate File was not Read Correctly\n");
        fatal();
    }
}
