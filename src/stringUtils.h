//////////////////////////////////////////////////////////
//                                                      //
//  stringUtils.h  --  Utility for string manipulation  //
//                                                      //
//////////////////////////////////////////////////////////


#pragma once
#include <string>
#include <vector>

namespace stringUtils
{
std::vector<std::string> split(std::string line);
bool checkIsDouble(std::string inputString, real& result);
}
