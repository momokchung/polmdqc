////////////////////////////////////////////////
//                                            //
//  readjson.h  --  read json to std::vector  //
//                                            //
////////////////////////////////////////////////


#pragma once
#include <string>
#include <vector>

std::vector<double> readVectorFromJson(const std::string& filename);
std::vector<std::vector<double>> readMatrixFromJson(const std::string& filename);
