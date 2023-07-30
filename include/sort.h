////////////////////////////////////////////////
//                                            //
//  sort.h  --  sorts and removes duplicates  //
//                                            //
////////////////////////////////////////////////


#pragma once
#include <vector>

template <typename T>
void sort(std::vector<T>& vector);

template <typename T>
void sortkey(std::vector<T>& vector, std::vector<int>& key);
