// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include <string>
#include <vector>

namespace stringUtils
{
////////////////////////////////////////////////////////
//                                                    //
//  stringUtils  --  Utility for string manipulation  //
//                                                    //
////////////////////////////////////////////////////////

std::vector<std::string> split(std::string line);
bool checkIsDouble(std::string inputString, real& result);
}
