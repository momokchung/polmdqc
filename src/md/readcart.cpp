// Author: Moses KJ Chung
// Year:   2023

#include "output.h"
#include "readcart.h"
#include "readxyz.h"

namespace polmdqc
{
////////////////////////////////////////////////////
//                                                //
//  readcart  --  input of Cartesian coordinates  //
//                                                //
////////////////////////////////////////////////////

// "readcart" gets a set of Cartesian coordinates from either
// a formatted or binary disk file

void readcart(std::ifstream& ffile, bool& first)
{
    // get next coordinates set from formatted or binary file
    if (archive) {
        readxyz(ffile);
    }
    else if (binary) {
        // readdcd(ffile,first);
    }
}
}
