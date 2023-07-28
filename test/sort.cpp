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
}
