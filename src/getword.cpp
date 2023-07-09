/////////////////////////////////////////////////////////
//                                                     //
//  getword.cpp  --  extract first word from a string  //
//                                                     //
/////////////////////////////////////////////////////////

// "getword" searches an input string for the first alphabetic
// character (A-Z or a-z); the region from this first character
// to the first blank space or separator is returned as a "word";
// if the actual word is too long, only the first part is returned

// variables and parameters:

// string    input character string to be searched
// word      output with the first word in the string
// next      input with first position of search string;
            // output with the position following word


#include "gettext.h"
#include "getword.h"
#include <iostream>
#include <sstream>

void getword(std::string& string, std::string& word, int& next)
{
    // get starting index
    int start = -1;
    int stringLength = string.length();
    for (int i = next; i < stringLength; i++) {
        char c = string[i];
        if (std::isalpha(c)) {
            start = i;
            break;
        }
    }

    // get ending index
    int end = -1;
    for (int i = start; i < stringLength; i++) {
        char c = string[i];
        if (c == ' ' or c == '\t' or c == ',' or c == ':' or c == ';') {
            end = i;
            break;
        }
    }

    if (start == -1 or end == -1) {
        word = "";
        return;
    }
    
    word = string.substr(start, end-start);

    next = end;
}
