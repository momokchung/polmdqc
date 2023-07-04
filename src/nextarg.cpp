////////////////////////////////////////////////////////
//                                                    //
//  nextarg.cpp  --  find next command line argument  //
//                                                    //
////////////////////////////////////////////////////////

// "nextarg" finds the next unused command line argument
// and returns it in the input character string


#include "argue.h"
#include "nextarg.h"

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
