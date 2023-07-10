#include "catch.hpp"
#include "sort.h"

TEST_CASE("sort") {
    std::vector<int> numbers = {5, 2, 8, 2, 1, 5, 9, 8, 2};
    std::vector<int> sortedNumbers = {1, 2, 5, 8, 9};
    sort<int>(numbers);
    REQUIRE(numbers == sortedNumbers);
}
