// Author: Moses KJ Chung
// Year:   2023

#include "upcase.h"

namespace polmdqc
{
//////////////////////////////////////////////////////
//                                                  //
//  upcase  --  uppercase all characters in string  //
//                                                  //
//////////////////////////////////////////////////////

void upcase(std::string& str)
{
    for (char& c : str) {
        c = std::toupper(c);
    }
}
}
