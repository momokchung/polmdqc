///////////////////////////////////////////////
//                                           //
//  ioUtils.h  --  Utility for input/output  //
//                                           //
///////////////////////////////////////////////


#pragma once
#include <string>
#include <vector>

namespace ioUtils{
void fileExists(std::string fileName);
int lineNumbers(std::string fileName);
void readlines(std::string fileName, std::vector<std::string>& lines);
}
