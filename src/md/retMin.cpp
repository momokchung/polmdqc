// Author: Moses KJ Chung
// Year:   2024

#include "retmin.h"

namespace polmdqc
{
////////////////////////////////////////////////////////
//                                                    //
//  retmin  --  return the minimum value of an array  //
//                                                    //
////////////////////////////////////////////////////////

// "retmin" return the minimum value of the array

template<typename T>
T retmin(MDQCArray<T>& arr)
{
    int size = arr.size();
    T min = arr[0];

    for (int i = 1; i < size; i++) {
        if (arr[i] < min) min = arr[i];
    }

    return min;
}

template real retmin(MDQCArray<real>& arr);
template int retmin(MDQCArray<int>& arr);
}
