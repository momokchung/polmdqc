#include "catch.hpp"
#include "precision.h"
#include "sort.h"
#include <set>

namespace polmdqc
{
TEST_CASE("sort") {
    int n;
    std::vector<int> numbers;
    std::vector<int> sortedNumbers;

    n = 9;
    numbers = {5, 2, 8, 2, 1, 5, 9, 8, 2};
    sortedNumbers = {1, 2, 5, 8, 9, 5, 9, 8, 2};
    sortUnique(n, numbers);
    REQUIRE(numbers == sortedNumbers);
    REQUIRE(n == 5);

    n = 5;
    numbers = {5, 2, 8, 2, 1, 5, 9, 8, 2};
    sortedNumbers = {5, 1, 2, 5, 8, 5, 9, 8, 2};
    sortUnique(n, numbers, 1);
    REQUIRE(numbers == sortedNumbers);
    REQUIRE(n == 4);

    n = 10;
    std::vector<real> doubles = {1.1, 2.2, -3.3, -2.2, 2.2, 3.3, -1.1, 0., 2.2, -3.3};
    std::vector<real> sortedDoubles = {-3.3, -2.2, -1.1, 0., 1.1, 2.2, 3.3, 0., 2.2, -3.3};
    sortUnique(n, doubles);
    REQUIRE(doubles == sortedDoubles);
    REQUIRE(n == 7);

    n = 7;
    doubles = {1.1, 2.2, -3.3, -2.2, 2.2, 3.3, -1.1, 0., 2.2, -3.3};
    sortedDoubles = {1.1, 2.2, -3.3, -2.2, -1.1, 0., 2.2, 3.3, 2.2, -3.3};
    sortUnique(n, doubles, 2);
    REQUIRE(doubles == sortedDoubles);
    REQUIRE(n == 6);

    numbers = {5, 2, 8, 2, 1, 5, 9, 8, 2};
    sortedNumbers = {1, 2, 2, 2, 5, 5, 8, 8, 9};
    std::vector<int> key;
    n = numbers.size();
    key.resize(n);
    sortKey(n, numbers, key);
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

    numbers = {5, 2, 8, 2, 1, 5, 9, 8, 2};
    sortedNumbers = {2, 2, 5, 8, 1, 5, 9, 8, 2};
    n = numbers.size();
    key.resize(n);
    sortKey(4, numbers, key);
    std::set<int> uniqueInts2(key.begin(), key.begin()+4);
    int numUniqueInts2 = uniqueInts2.size();
    REQUIRE(numbers == sortedNumbers);
    REQUIRE(numUniqueInts2 == 4);
    REQUIRE((key[0] == 1 or key[0] == 3));
    REQUIRE((key[1] == 1 or key[1] == 3));
    REQUIRE(key[2] == 0);
    REQUIRE(key[3] == 2);

    doubles = {1.1, 2.2, -3.3, -2.2, 2.2, 3.3, -1.1, 0., 2.2, -3.3};
    sortedDoubles = {-3.3, -3.3, -2.2, -1.1, 0., 1.1, 2.2, 2.2, 2.2, 3.3};
    n = doubles.size();
    key.resize(n);
    sortKey(n, doubles, key);
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

    doubles = {1.1, 2.2, -3.3, -2.2, 2.2, 3.3, -1.1, 0., 2.2, -3.3};
    sortedDoubles = {-3.3, -2.2, 1.1, 2.2, 2.2, 3.3, -1.1, 0., 2.2, -3.3};
    n = doubles.size();
    key.resize(n);
    sortKey(6, doubles, key);
    std::set<int> uniqueDoubles2(key.begin(), key.begin()+6);
    int numUniqueDoubles2 = uniqueDoubles2.size();
    REQUIRE(doubles == sortedDoubles);
    REQUIRE(numUniqueDoubles2 == 6);
    REQUIRE(key[0] == 2);
    REQUIRE(key[1] == 3);
    REQUIRE(key[2] == 0);
    REQUIRE((key[3] == 1 or key[3] == 4));
    REQUIRE((key[4] == 1 or key[4] == 4));
    REQUIRE(key[5] == 5);
}
}
