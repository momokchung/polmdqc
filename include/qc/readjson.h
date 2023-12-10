// Author: Moses KJ Chung
// Year:   2023

#pragma once
#include <string>
#include <vector>

//////////////////////////////////////////////
//                                          //
//  readjson  --  read json to std::vector  //
//                                          //
//////////////////////////////////////////////

std::vector<double> readVectorFromJson(const std::string& filename);
std::vector<std::vector<double>> readMatrixFromJson(const std::string& filename);
