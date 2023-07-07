////////////////////////////////////////////////////////////
//                                                        //
//  stringUtils.cpp  --  Utility for string manipulation  //
//                                                        //
////////////////////////////////////////////////////////////


#include "init.h"
#include <sstream>
#include <stdlib.h>
#include <string>
#include <vector>

namespace stringUtils
{
/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//  vector<string> split  --  splits string line by into a vector<string>  //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

std::vector<std::string> split(std::string line)
{
    std::vector<std::string> words{};
    std::string word;
    std::istringstream iss(line);
    while (iss >> word)
    {
        words.push_back(word);
    }
    return words;
}


//////////////////////////////////////////////////////////////////////////////////////
//                                                                                  //
//  bool checkIsDouble  --  checks if string is double, if true converts to double  //
//                                                                                  //
//////////////////////////////////////////////////////////////////////////////////////

bool checkIsDouble(std::string inputString, real& result) {
    char* end;
    result = std::strtod(inputString.c_str(), &end);
    if (end == inputString.c_str() || *end != '\0') return false;
    return true;
}
}
