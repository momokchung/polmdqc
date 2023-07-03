///////////////////////////////////////////////////////
//                                                   //
//  command.cpp  --  get any command line arguments  //
//                                                   //
///////////////////////////////////////////////////////

// "command" uses the standard Unix-like iargc/getarg routines
// to get the number and values of arguments specified on the
// command line at program runtime


#include "command.h"
#include "argue.h"
#include <iostream>
#include <string>
#include <cstring>
#include <cctype>

std::string arg[maxarg];
int narg;
bool listarg[maxarg];

void command(int argc, char** argv)
{
    // initialize command line arguments as blank strings
    std::string blank = "                    ";
    for (int i = 0; i < maxarg; ++i) {
        arg[i] = blank + blank + blank;
    }

    // get the number of arguments and store each in a string
    narg = std::min(argc, maxarg);
    for (int i = 0; i < narg; ++i) {
        arg[i] = argv[i];
    }

    // mark the command line options as unuseable for input
    listarg[0] = false;
    for (int i = 1; i < narg; ++i) {
        listarg[i] = true;
    }
    for (int i = 1; i <= narg; ++i) {
        char letter = arg[i][0];
        if (letter == '-') {
            letter = arg[i][1];
            letter = std::toupper(static_cast<unsigned char>(letter));
            if (letter >= 'A' && letter <= 'Z') {
                listarg[i] = false;
                listarg[i+1] = false;
            }
        }
    }
}
