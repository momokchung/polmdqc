// Author: Moses KJ Chung
// Year:   2024

#include "retmax.h"

namespace polmdqc
{
////////////////////////////////////////////////////////
//                                                    //
//  retmax  --  return the maximum value of an array  //
//                                                    //
////////////////////////////////////////////////////////

// "retmax" return the maximum value of the array

template<typename T>
T retmax(MDQCArray<T>& arr)
{
    int size = arr.size();
    T max = arr[0];

    for (int i = 1; i < size; i++) {
        if (arr[i] > max) max = arr[i];
    }

    return max;
}

template real retmax(MDQCArray<real>& arr);
template int retmax(MDQCArray<int>& arr);
}
