// Author: Moses KJ Chung
// Year:   2023

#include "getline.h"
#include <iostream>

namespace polmdqc
{
///////////////////////////////////////////////
//                                           //
//  gettext  --  extract text from a string  //
//                                           //
///////////////////////////////////////////////

// "getline" trims whitespace of a given string

void getline(std::string& string)
{
    std::size_t start = string.find_first_not_of(" \t\n\r");
    std::size_t end = string.find_last_not_of(" \t\n\r");
  
    if (start == std::string::npos || end == std::string::npos) {
        string = "";
        return;
    }

    string = string.substr(start, end - start + 1);
}
}
