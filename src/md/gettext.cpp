// Author: Moses KJ Chung
// Year:   2023

#include "gettext.h"
#include <iostream>
#include <sstream>

namespace polmdqc
{
///////////////////////////////////////////////
//                                           //
//  gettext  --  extract text from a string  //
//                                           //
///////////////////////////////////////////////

// "gettext" searches an input string for the first string of
// non-blank characters; the region from a non-blank character
// to the first space or tab is returned as "text"; if the
// actual text is too long, only the first part is returned

// variables and parameters:

// string    input character string to be searched
// text      output with the first text string found
// next      input with first position of search string;
//              output with the position following text

void gettext(std::string& string, std::string& text, int& next)
{
    std::string stringSub = string.substr(next);
    std::istringstream iss(stringSub);
    std::string firstWord = "";
    iss >> firstWord;
    text = firstWord;
    int numSpaces = 0;
    for (char ch : stringSub) {
        if (ch == ' ' or ch == '\t')
            numSpaces++;
        else
            break;
    }
    next += text.length() + numSpaces;
}
}
