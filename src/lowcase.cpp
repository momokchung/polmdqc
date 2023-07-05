///////////////////////////////////////////////////////////
//                                                       //
//  lowcase.cpp  --  lowercase all characters in string  //
//                                                       //
///////////////////////////////////////////////////////////


#include "lowcase.h"

void lowcase(std::string& str)
{
    for (char& c : str) {
        c = std::tolower(c);
    }
}
