// Author: Moses KJ Chung
// Year:   2023

#include "getstring.h"
#include <iostream>

namespace polmdqc
{
///////////////////////////////////////////////////
//                                               //
//  getstring  --  extract double quoted string  //
//                                               //
///////////////////////////////////////////////////

// "getstring" searches for a quoted text string within an input
// character string; the region between the first and second
// double quote is returned as the "text"; if the actual text is
// too long, only the first part is returned
//
// variables and parameters:
//
// string    input character string to be searched
// text      the quoted text found in the input string
// next      input with first position of search string;
//             output with the position following text

void getstring(std::string& string, std::string& text, int& next)
{
    std::string newString = string.substr(next);
    size_t firstQuotePos = newString.find("\"");
    size_t secondQuotePos = newString.find("\"", firstQuotePos + 1);

    if (firstQuotePos != std::string::npos and secondQuotePos != std::string::npos) {
        text = newString.substr(firstQuotePos + 1, secondQuotePos - firstQuotePos - 1);
        next += secondQuotePos + 1;
    }
    else {
        text = "";
    }
    
}
}
