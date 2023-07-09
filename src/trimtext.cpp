///////////////////////////////////////////////////////
//                                                   //
//  trimtext.cpp  --  find last non-blank character  //
//                                                   //
///////////////////////////////////////////////////////

// "trimtext" finds and returns the location of the last
// non-blank character before the first null character in
// an input text string; the function returns -1 if no
// such character is found


#include "trimtext.h"

int trimtext(std::string& string)
{
    int lastIndex = -1;
    for (int i = string.length()-1; i >= 0; i--) {
        if (!std::isspace(string[i])) {
            lastIndex = i;
            break;
        }
    }
    return lastIndex;
}
