// Author: Moses KJ Chung
// Year:   2023

#include "argue.h"
#include "nextarg.h"

namespace polmdqc
{
////////////////////////////////////////////////////
//                                                //
//  nextarg  --  find next command line argument  //
//                                                //
////////////////////////////////////////////////////

// "nextarg" finds the next unused command line argument
// and returns it in the input character string

void nextarg(std::string& string, bool& exist)
{
    // initialize the command argument as a blank string
    string = "          ";
    exist = false;

    // get the next command line argument and mark it as used
    if (narg != 0) {
        for (int i = 1; i <= narg; i++) {
            if (listarg[i]) {
                listarg[i] = false;
                string = arg[i];
                exist = true;
                break;
            }
        }
    }
}
}
