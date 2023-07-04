/////////////////////////////////////////////////////////
//                                                     //
//  getcart.cpp  --  get a Cartesian coordinates file  //
//                                                     //
/////////////////////////////////////////////////////////

// "getcart" asks for a Cartesian coordinate file name, then
// reads the formatted or binary coordinates file


#include "basefile.h"
#include "getcart.h"
#include "nextarg.h"
#include <string>

void getcart(int &ixyz)
{
    std::string xyzfile;
    bool exist;
    // try to get a filename from the command line arguments
    nextarg (xyzfile,exist);
    if (exist) {
        basefile(xyzfile);
        // call suffix (xyzfile,'xyz','old')
        // inquire (file=xyzfile,exist=exist)
    }

}
