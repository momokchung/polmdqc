// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include <fstream>

namespace polmdqc
{
////////////////////////////////////////////////////
//                                                //
//  readcart  --  input of Cartesian coordinates  //
//                                                //
////////////////////////////////////////////////////

void readcart(std::ifstream& ffile, bool& first);
}
