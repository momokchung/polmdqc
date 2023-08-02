#include "catch.hpp"
#include "sort.h"
#include <set>

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
    sortedNumbers = {5, 1, 2, 5, 8, 9, 2};
    doubles = {1.1, 2.2, -3.3, -2.2, 2.2, 3.3, -1.1, 0., 2.2, -3.3};
    sortedDoubles = {1.1, -3.3, -2.2, -1.1, 2.2, 3.3, 0., 2.2, -3.3};
    sort(numbers, 1, 8);
    REQUIRE(numbers == sortedNumbers);
    sort(doubles, 1, 7);
    REQUIRE(doubles == sortedDoubles);
    numbers = {5, 2, 8, 2, 1, 5, 9, 8, 2};
    sortedNumbers = {1, 2, 2, 2, 5, 5, 8, 8, 9};
    std::vector<int> key;
    sortkey(numbers, key);
    std::set<int> uniqueInts(key.begin(), key.end());
    int numUniqueInts = uniqueInts.size();
    REQUIRE(numbers == sortedNumbers);
    REQUIRE(numUniqueInts == 9);
    REQUIRE(key[0] == 4);
    REQUIRE((key[1] == 1 or key[1] == 3 or key[1] == 8));
    REQUIRE((key[2] == 1 or key[2] == 3 or key[2] == 8));
    REQUIRE((key[3] == 1 or key[3] == 3 or key[3] == 8));
    REQUIRE((key[4] == 0 or key[4] == 5));
    REQUIRE((key[5] == 0 or key[5] == 5));
    REQUIRE((key[6] == 2 or key[6] == 7));
    REQUIRE((key[7] == 2 or key[7] == 7));
    REQUIRE(key[8] == 6);
    doubles = {1.1, 2.2, -3.3, -2.2, 2.2, 3.3, -1.1, 0., 2.2, -3.3};
    sortedDoubles = {-3.3, -3.3, -2.2, -1.1, 0., 1.1, 2.2, 2.2, 2.2, 3.3};
    sortkey(doubles, key);
    std::set<int> uniqueDoubles(key.begin(), key.end());
    int numUniqueDoubles = uniqueDoubles.size();
    REQUIRE(doubles == sortedDoubles);
    REQUIRE(numUniqueDoubles == 10);
    REQUIRE((key[0] == 2 or key[0] == 9));
    REQUIRE((key[1] == 2 or key[1] == 9));
    REQUIRE(key[2] == 3);
    REQUIRE(key[3] == 6);
    REQUIRE(key[4] == 7);
    REQUIRE(key[5] == 0);
    REQUIRE((key[6] == 1 or key[6] == 4 or key[6] == 8));
    REQUIRE((key[7] == 1 or key[7] == 4 or key[7] == 8));
    REQUIRE((key[8] == 1 or key[8] == 4 or key[8] == 8));
    REQUIRE(key[9] == 5);
}
