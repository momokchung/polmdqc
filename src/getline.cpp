///////////////////////////////////////////////////
//                                               //
//  gettext.cpp  --  extract text from a string  //
//                                               //
///////////////////////////////////////////////////

// "getline" trims whitespace of a given string


#include "getline.h"
#include <iostream>

std::string getline(std::string& string)
{
    std::size_t start = string.find_first_not_of(" \t\n\r");
    std::size_t end = string.find_last_not_of(" \t\n\r");
  
    if (start == std::string::npos || end == std::string::npos)
        return "";
  
    return string.substr(start, end - start + 1);
}
