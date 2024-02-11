// Author: Moses KJ Chung
// Year:   2024

#include "rootdir.h"
#include <filesystem>

namespace polmdqc
{
////////////////////////////////////////////////////
//                                                //
//  rootdir  --  get file path of root directory  //
//                                                //
////////////////////////////////////////////////////

// "rootdir" returns the path of polmdqc root directory

std::string rootdir()
{
    std::filesystem::path filePath = __FILE__;
    std::filesystem::path directoryPath = filePath.parent_path().parent_path();
    // std::string path = directoryPath;
    return directoryPath;
}
}
