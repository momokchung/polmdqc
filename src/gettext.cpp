///////////////////////////////////////////////////
//                                               //
//  gettext.cpp  --  extract text from a string  //
//                                               //
///////////////////////////////////////////////////

// "gettext" searches an input string for the first string of
// non-blank characters; the region from a non-blank character
// to the first space or tab is returned as "text"; if the
// actual text is too long, only the first part is returned

// variables and parameters:

// string    input character string to be searched
// text      output with the first text string found
// next      input with first position of search string;
//              output with the position following text


#include "gettext.h"
#include <iostream>
#include <sstream>

void gettext(std::string& string, std::string& text, int& next)
{
    std::istringstream iss(string.substr(next));
    std::string firstWord;
    iss >> firstWord;
    text = firstWord;
    int numSpaces = 0;
    for (char ch : string) {
        if (ch == ' ')
            numSpaces++;
        else
            break;
    }
    next = text.length() + numSpaces;
}
