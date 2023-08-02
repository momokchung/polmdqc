//////////////////////////////////////////////////
//                                              //
//  sort.cpp  --  sorts and removes duplicates  //
//                                              //
//////////////////////////////////////////////////

// "sortUnique" sorts and removes duplicates from a std::vector object
// and returns the number of unique objects


#include "sort.h"
#include <algorithm>

template <typename T>
void sortUnique(int& n, std::vector<T>& vector, size_t startIndex, size_t endIndex)
{
    // Invalid subset range, do nothing
    if (startIndex >= endIndex || startIndex >= vector.size()) {
        return;
    }

    // Adjust endIndex if it exceeds the vector size
    endIndex = std::min(endIndex, vector.size());

    // Sort the subset within the given range
    std::sort(vector.begin() + startIndex, vector.begin() + endIndex);

    // remove duplicates
    auto last = std::unique(vector.begin() + startIndex, vector.begin() + endIndex);
    vector.erase(last, vector.begin() + endIndex);
    n = vector.size();
}

template void sortUnique<int>(int& n, std::vector<int>& vector, size_t startIndex, size_t endIndex);
template void sortUnique<double>(int& n, std::vector<double>& vector, size_t startIndex, size_t endIndex);


// "sortKey" sorts and returns a key into the original ordering

template <typename T>
void sortKey(std::vector<T>& vector, std::vector<int>& key)
{
    int n = vector.size();
    key.resize(n);

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

template void sortKey<int>(std::vector<int>& vector, std::vector<int>& key);
template void sortKey<double>(std::vector<double>& vector, std::vector<int>& key);
