////////////////////////////////////////////////
//                                            //
//  sort.h  --  sorts and removes duplicates  //
//                                            //
////////////////////////////////////////////////


#pragma once
#include <vector>

template <typename T>
void sortUnique(int& n, std::vector<T>& vector, int startIndex=0);

template <typename T>
void sortKey(size_t n, std::vector<T>& vector, std::vector<int>& key);
