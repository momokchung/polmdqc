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
    gettext(string, word, next);

    // get starting index
    int start = 0;
    for (char c : word) {
        if (std::isalpha(c)) {
            break;
        }
        start++;
    }

    // get ending index
    int end = 0;
    for (char c : word.substr(start)) {
        if (!std::isalnum(c)) {
            break;
        }
        end++;
    }

    word = word.substr(start, end);
}
