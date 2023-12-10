// Author: Moses KJ Chung
// Year:   2023

#include "inquire.h"
#include <fstream>

namespace polmdqc
{
/////////////////////////////////////////
//                                     //
//  inquire  --  check if file exists  //
//                                     //
/////////////////////////////////////////

bool inquireFile(std::string& str)
{
    std::ifstream f(str);
    return f.good();
}

bool inquireUnit(std::ifstream& unit)
{
    return unit.is_open();
}
}
