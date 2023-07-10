//////////////////////////////////////////////////
//                                              //
//  sort.cpp  --  sorts and removes duplicates  //
//                                              //
//////////////////////////////////////////////////

// "sort" sorts and removes duplicates from a std::vector object


#include "sort.h"
#include <algorithm>

template <typename T>
void sort(std::vector<T>& vector)
{
    // remove duplicates
    std::sort(vector.begin(), vector.end());
    auto last = std::unique(vector.begin(), vector.end());
    vector.erase(last, vector.end());

    // // sort the vector
    // std::sort(vector.begin(), vector.end());
}

template void sort<int>(std::vector<int>& vector);
template void sort<double>(std::vector<double>& vector);
