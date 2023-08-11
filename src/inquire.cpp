/////////////////////////////////////////////
//                                         //
//  inquire.cpp  --  check if file exists  //
//                                         //
/////////////////////////////////////////////


#include "inquire.h"
#include <fstream>

bool inquireFile(std::string& str)
{
    std::ifstream f(str);
    return f.good();
}

bool inquireUnit(std::ifstream& unit)
{
    return unit.is_open();
}
