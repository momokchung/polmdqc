//////////////////////////////////////////////////////////
//                                                      //
//  upcase.cpp  --  uppercase all characters in string  //
//                                                      //
//////////////////////////////////////////////////////////


#include "upcase.h"

void upcase(std::string& str)
{
    for (char& c : str) {
        c = std::toupper(c);
    }
}
