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
}

template void sort<int>(std::vector<int>& vector);
template void sort<double>(std::vector<double>& vector);


// "sortkey" sorts and returns a key into the original ordering

template <typename T>
void sortkey(std::vector<T>& vector, std::vector<int>& key)
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

template void sortkey<int>(std::vector<int>& vector, std::vector<int>& key);
template void sortkey<double>(std::vector<double>& vector, std::vector<int>& key);
