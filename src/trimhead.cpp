/////////////////////////////////////////////////////////
//                                                     //
//  trimhead.cpp  --  remove spaces before first text  //
//                                                     //
/////////////////////////////////////////////////////////

// "trimhead" removes blank spaces before the first non-blank
// character in a text string by shifting the string to the left


#include "trimhead.h"

void trimhead(std::string& string)
{
    size_t firstNonBlank = string.find_first_not_of(' ');

    if (firstNonBlank != std::string::npos && firstNonBlank > 0) {
        string = string.substr(firstNonBlank);
    }
}
