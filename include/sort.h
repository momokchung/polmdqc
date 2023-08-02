////////////////////////////////////////////////
//                                            //
//  sort.h  --  sorts and removes duplicates  //
//                                            //
////////////////////////////////////////////////


#pragma once
#include <limits>
#include <vector>

template <typename T>
void sortUnique(int& n, std::vector<T>& vector, size_t startIndex=0, size_t endIndex=std::numeric_limits<size_t>::max());

template <typename T>
void sortkey(std::vector<T>& vector, std::vector<int>& key);
