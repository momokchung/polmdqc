// Author: Moses KJ Chung
// Year:   2023

#include "getnumb.h"
#include <sstream>

namespace polmdqc
{
/////////////////////////////////////////////////////
//                                                 //
//  getnumb  --  extract an integer from a string  //
//                                                 //
/////////////////////////////////////////////////////

// "getnumb" searches an input string from left to right for an
// integer and puts the numeric value in "number"; returns zero
// with "next" unchanged if no integer value is found
//
// variables and parameters:
//
// string    input character string to be searched
// number    output with the first integer in the string
// next      input with first position of search string;
//             output with the position following the number

void getnumb(std::string& string, int& number, int& next)
{
    number = 0;
    std::string newString = string.substr(next);
    std::istringstream iss(newString);

    if (!(iss >> number)) return;

    std::streampos pos = iss.tellg();
    int posInt = static_cast<int>(pos);
    int length = newString.length();
    int addLength = (posInt == -1 ? length : posInt);

    next += addLength;
}
}
