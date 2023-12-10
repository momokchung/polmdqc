// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include <fstream>
#include <string>

namespace polmdqc
{
/////////////////////////////////////////
//                                     //
//  inquire  --  check if file exists  //
//                                     //
/////////////////////////////////////////

bool inquireFile(std::string& str);

bool inquireUnit(std::ifstream& unit);
}
