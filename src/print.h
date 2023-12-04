//////////////////////////////////////////////
//                                          //
//  print.h  --  various printing routines  //
//                                          //
//////////////////////////////////////////////


#pragma once
#include "init.h"
#include <string>
#include <vector>

namespace print
{
void printMatrix(const std::vector<std::vector<real>>& matrix, std::string header);
void printVMatrix(const std::vector<real>& vmatrix, int rowN, int colN, std::string header);
}
