//////////////////////////////////////////////////
//                                              //
//  sort.cpp  --  sorts and removes duplicates  //
//                                              //
//////////////////////////////////////////////////

// "sortUnique" sorts and removes duplicates from a std::vector
// object and returns the number of unique objects


#include "sort.h"
#include <algorithm>
#include <unordered_set>

template <typename T>
void sortUnique(int& n, std::vector<T>& vector, int startIndex)
{
    std::unordered_set<T> s;
    for (int i = startIndex; i < startIndex + n; i++) {
        s.insert(vector[i]);
    }
    int setSize = s.size();

    size_t start = startIndex;
    for (auto it = s.begin(); it != s.end(); ++it) {
        vector[start] = *it;
        ++start;
    }

    std::sort(vector.begin()+startIndex, vector.begin()+startIndex+setSize);
    n = setSize;
}

template void sortUnique<int>(int& n, std::vector<int>& vector, int startIndex);
template void sortUnique<double>(int& n, std::vector<double>& vector, int startIndex);


// "sortKey" sorts and returns a key into the original ordering

template <typename T>
void sortKey(size_t n, std::vector<T>& vector, std::vector<int>& key)
{
    // Invalid range, do nothing
    if (n > vector.size()) return;

    // Create a temporary vector with pairs of (value, index)
    std::vector<std::pair<T,int>> temp;
    for (int i = 0; i < n; ++i) {
        temp.push_back({vector[i], i});
    }

    // Sort the temporary vector based on the values
    std::sort(temp.begin(), temp.end(), [](const std::pair<T,int>& a, const std::pair<T,int>& b) {
        return a.first < b.first;
    });

    // Update the key vector with the original indices in the sorted order
    key.clear();
    for (const auto& pair : temp) {
        key.push_back(pair.second);
    }

    // Update the original vector with the sorted values
    for (int i = 0; i < n; ++i) {
        vector[i] = temp[i].first;
    }
}

template void sortKey<int>(size_t n, std::vector<int>& vector, std::vector<int>& key);
template void sortKey<double>(size_t n, std::vector<double>& vector, std::vector<int>& key);
