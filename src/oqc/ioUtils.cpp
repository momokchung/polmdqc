// Author: Moses KJ Chung
// Year:   2023

#include "ioUtils.h"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace ioUtils{
////////////////////////////////////////////////////////////
//                                                        //
//  fileExists  --  exits program if file does not exist  //
//                                                        //
////////////////////////////////////////////////////////////

void fileExists(std::string fileName)
{
    std::ifstream f(fileName);
    if (!f)
    {
        std::cerr << "Unable to open " << fileName << "\n";
        std::exit(1);
    }
    f.close();
}

///////////////////////////////////////////////////////////////
//                                                           //
//  lineNumbers  --  returns total number of nonempty lines  //
//                                                           //
///////////////////////////////////////////////////////////////

int lineNumbers(std::string fileName)
{
    std::ifstream f(fileName);
    std::string line;
    int lineN = 0;
    while (std::getline(f, line))
    {
        if (line.length() != 0)
        {
            lineN += 1;
        }
    }
    f.close();
    return lineN;
}

/////////////////////////////////////////////////////////////////
//                                                             //
//  readlines  --  returns nonempty lines as a vector<string>  //
//                                                             //
/////////////////////////////////////////////////////////////////

void readlines(std::string fileName, std::vector<std::string>& lines)
{
    std::ifstream f(fileName);
    std::string line;
    int i = 0;
    while (std::getline(f, line))
    {
        if (line.length() != 0)
        {
            std::transform(line.begin(), line.end(), line.begin(), toupper);
            lines[i] = line;
            i += 1;
        }
    }
    f.close();
}
}
