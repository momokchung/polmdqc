// Author: Moses KJ Chung
// Year:   2023

#include "init.h"
#include <iostream>
#include <filesystem>

namespace init
{
////////////////////////////
//                        //
//  init  --  initialize  //
//                        //
////////////////////////////

std::string cwd;

namespace fs = std::filesystem;
void init(char** argv)
{
    try {
        // Replace "your_file.txt" with the actual file name or path
        fs::path filePath = __FILE__;

        // Get the directory containing the file
        fs::path directoryPath = filePath.parent_path().parent_path();

        cwd = directoryPath;
    } catch (const fs::filesystem_error& ex) {
        std::cerr << "Filesystem error: " << ex.what() << std::endl;
    }
}
}