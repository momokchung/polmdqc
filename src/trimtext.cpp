///////////////////////////////////////////////////////
//                                                   //
//  trimtext.cpp  --  find last non-blank character  //
//                                                   //
///////////////////////////////////////////////////////

// "trimtext" finds and returns the location of the last
// non-blank character before the first null character in
// an input text string; the function returns zero if no
// such character is found


#include "trimtext.h"

int trimtext(std::string& string)
{
    // move forward through the string, one character
    // at a time, looking for first null character
    int trimtext = 0;
    int size = string.length();
    char null = '\0';
    int last = size;
    for (int i = 0; i < size; i++) {
        if (string[i] == null) {
            last = i;
            break;
        }
    }

    // move backward through the string, one character
    // at a time, looking for first non-blank character
    for (int i = last; i >= 0; i--) {
        if (string[i] > ' ') {
            trimtext = i;
            break;
        }
    }

    return trimtext + 1;
}
