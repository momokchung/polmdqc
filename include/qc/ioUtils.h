// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include <string>
#include <vector>

namespace ioUtils{
/////////////////////////////////////////////
//                                         //
//  ioUtils  --  Utility for input/output  //
//                                         //
/////////////////////////////////////////////

void fileExists(std::string fileName);
int lineNumbers(std::string fileName);
void readlines(std::string fileName, std::vector<std::string>& lines);
}
