// Author: Moses KJ Chung
// Year:   2023

#include "getline.h"
#include "justify.h"

namespace polmdqc
{
//////////////////////////////////////////////////////
//                                                  //
//  justify  --  convert string to right justified  //
//                                                  //
//////////////////////////////////////////////////////

// "justify" converts a text string to right justified format
// with leading blank spaces

void justify(std::string& string, int totalLength)
{
    getline(string);
    int stringLength = string.length();
    if (totalLength <= stringLength)
    {
        return;
    }

    int numSpaces = totalLength - stringLength;
    std::string spaces(numSpaces, ' ');
    
    string = spaces + string;
}
}
