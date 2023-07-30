#include "catch.hpp"
#include "sort.h"

TEST_CASE("sort") {
    std::vector<int> numbers = {5, 2, 8, 2, 1, 5, 9, 8, 2};
    std::vector<int> sortedNumbers = {1, 2, 5, 8, 9};
    std::vector<double> doubles = {1.1, 2.2, -3.3, -2.2, 2.2, 3.3, -1.1, 0., 2.2, -3.3};
    std::vector<double> sortedDoubles = {-3.3, -2.2, -1.1, 0., 1.1, 2.2, 3.3};
    sort(numbers);
    REQUIRE(numbers == sortedNumbers);
    sort(doubles);
    REQUIRE(doubles == sortedDoubles);
    numbers = {5, 2, 8, 2, 1, 5, 9, 8, 2};
    sortedNumbers = {1, 2, 2, 2, 5, 5, 8, 8, 9};
    std::vector<int> key;
    std::vector<int> sortedKey = {4, 1, 8, 3, 5, 0, 7, 2, 6};
    sortkey(numbers, key);
    REQUIRE(numbers == sortedNumbers);
    REQUIRE(key == sortedKey);
    doubles = {1.1, 2.2, -3.3, -2.2, 2.2, 3.3, -1.1, 0., 2.2, -3.3};
    sortedDoubles = {-3.3, -3.3, -2.2, -1.1, 0., 1.1, 2.2, 2.2, 2.2, 3.3};
    sortedKey = {9, 2, 3, 6, 7, 0, 4, 1, 8, 5};
    sortkey(doubles, key);
    REQUIRE(doubles == sortedDoubles);
    REQUIRE(key == sortedKey);
}
