//////////////////////////////////////////////////////////
//                                                      //
//  justify.cpp  --  convert string to right justified  //
//                                                      //
//////////////////////////////////////////////////////////

// "justify" converts a text string to right justified format
// with leading blank spaces


#include "getline.h"
#include "justify.h"

void justify(std::string& string, int totalLength)
{
    int stringLength = string.length();
    if (totalLength <= stringLength)
    {
        return;
    }

    int firstNonBlank = string.find_first_not_of(' ');
    int numSpaces = totalLength - stringLength + firstNonBlank;
    std::string spaces(numSpaces, ' ');
    
    string = spaces + string.substr(firstNonBlank);
}
