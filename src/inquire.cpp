/////////////////////////////////////////////
//                                         //
//  inquire.cpp  --  check if file exists  //
//                                         //
/////////////////////////////////////////////


#include "inquire.h"
#include <fstream>

bool inquire(std::string& str)
{
    std::ifstream f(str);
    return f.good();
}
