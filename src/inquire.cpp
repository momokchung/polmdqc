/////////////////////////////////////////////
//                                         //
//  inquire.cpp  --  check if file exists  //
//                                         //
/////////////////////////////////////////////


#include "inquire.h"
#include <fstream>

bool inquire(std::string& str)
{
    std::ifstream f(str.c_str());
    return f.good();
}
