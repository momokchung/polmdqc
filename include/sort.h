////////////////////////////////////////////////
//                                            //
//  sort.h  --  sorts and removes duplicates  //
//                                            //
////////////////////////////////////////////////


#pragma once
#include <vector>

template <typename T>
void sort(std::vector<T>& vector, size_t startIndex=0, size_t endIndex=std::numeric_limits<size_t>::max());

template <typename T>
void sortkey(std::vector<T>& vector, std::vector<int>& key);
