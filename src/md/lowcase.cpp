// Author: Moses KJ Chung
// Year:   2023

#include "lowcase.h"

namespace polmdqc
{
///////////////////////////////////////////////////////
//                                                   //
//  lowcase  --  lowercase all characters in string  //
//                                                   //
///////////////////////////////////////////////////////

void lowcase(std::string& str)
{
    for (char& c : str) {
        c = std::tolower(c);
    }
}
}
