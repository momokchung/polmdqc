#include "catch.hpp"
#include "precision.h"
#include "sort.h"
#include <set>

namespace polmdqc
{
TEST_CASE("sortUnique-1") {
    int n;
    std::vector<int> numbers;
    std::vector<int> sortedNumbers;
    std::vector<real> doubles;
    std::vector<real> sortedDoubles;

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
    doubles = {1.1, 2.2, -3.3, -2.2, 2.2, 3.3, -1.1, 0., 2.2, -3.3};
    sortedDoubles = {-3.3, -2.2, -1.1, 0., 1.1, 2.2, 3.3, 0., 2.2, -3.3};
    sortUnique(n, doubles);
    REQUIRE(doubles == sortedDoubles);
    REQUIRE(n == 7);

    n = 7;
    doubles = {1.1, 2.2, -3.3, -2.2, 2.2, 3.3, -1.1, 0., 2.2, -3.3};
    sortedDoubles = {1.1, 2.2, -3.3, -2.2, -1.1, 0., 2.2, 3.3, 2.2, -3.3};
    sortUnique(n, doubles, 2);
    REQUIRE(doubles == sortedDoubles);
    REQUIRE(n == 6);
}

TEST_CASE("sortUnique-2") {
    int n;
    int tn;
    std::vector<int> numbers;
    int* narr;
    std::vector<int> sortedNumbers;
    std::vector<real> doubles;
    real* darr;
    std::vector<real> sortedDoubles;

    n = 9;
    tn = 9;
    narr = new int[tn];
    numbers = {5, 2, 8, 2, 1, 5, 9, 8, 2};
    for (int i = 0; i < tn; i++) narr[i] = numbers[i];
    sortedNumbers = {1, 2, 5, 8, 9, 5, 9, 8, 2};
    sortUnique(n, narr);
    for (int i = 0; i < tn; i++) REQUIRE(narr[i] == sortedNumbers[i]);
    REQUIRE(n == 5);
    delete[] narr;

    n = 5;
    tn = 9;
    narr = new int[tn];
    numbers = {5, 2, 8, 2, 1, 5, 9, 8, 2};
    for (int i = 0; i < tn; i++) narr[i] = numbers[i];
    sortedNumbers = {5, 1, 2, 5, 8, 5, 9, 8, 2};
    sortUnique(n, narr, 1);
    for (int i = 0; i < tn; i++) REQUIRE(narr[i] == sortedNumbers[i]);
    REQUIRE(n == 4);
    delete[] narr;

    n = 10;
    tn = 10;
    darr = new real[tn];
    doubles = {1.1, 2.2, -3.3, -2.2, 2.2, 3.3, -1.1, 0., 2.2, -3.3};
    for (int i = 0; i < tn; i++) darr[i] = doubles[i];
    sortedDoubles = {-3.3, -2.2, -1.1, 0., 1.1, 2.2, 3.3, 0., 2.2, -3.3};
    sortUnique(n, darr);
    for (int i = 0; i < tn; i++) REQUIRE(darr[i] == sortedDoubles[i]);
    REQUIRE(n == 7);
    delete[] darr;

    n = 7;
    tn = 10;
    darr = new real[tn];
    doubles = {1.1, 2.2, -3.3, -2.2, 2.2, 3.3, -1.1, 0., 2.2, -3.3};
    for (int i = 0; i < tn; i++) darr[i] = doubles[i];
    sortedDoubles = {1.1, 2.2, -3.3, -2.2, -1.1, 0., 2.2, 3.3, 2.2, -3.3};
    sortUnique(n, darr, 2);
    for (int i = 0; i < tn; i++) REQUIRE(darr[i] == sortedDoubles[i]);
    REQUIRE(n == 6);
    delete[] darr;
}

TEST_CASE("sortKey-1") {
    int n;
    std::vector<int> key;
    std::vector<int> numbers;
    std::vector<int> sortedNumbers;
    std::vector<real> doubles;
    std::vector<real> sortedDoubles;

    numbers = {5, 2, 8, 2, 1, 5, 9, 8, 2};
    sortedNumbers = {1, 2, 2, 2, 5, 5, 8, 8, 9};
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

TEST_CASE("sortKey-2") {
    int n;
    int* key;
    std::vector<int> keyv;
    std::vector<int> numbers;
    std::vector<int> sortedNumbers;
    std::vector<real> doubles;
    std::vector<real> sortedDoubles;

    numbers = {5, 2, 8, 2, 1, 5, 9, 8, 2};
    sortedNumbers = {1, 2, 2, 2, 5, 5, 8, 8, 9};
    n = numbers.size();
    key = new int[n];
    keyv.resize(n);
    sortKey(n, numbers, key);
    for (int i = 0; i < n; i++) keyv[i] = key[i];
    std::set<int> uniqueInts(keyv.begin(), keyv.end());
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
    delete[] key;

    numbers = {5, 2, 8, 2, 1, 5, 9, 8, 2};
    sortedNumbers = {2, 2, 5, 8, 1, 5, 9, 8, 2};
    n = numbers.size();
    key = new int[n];
    keyv.resize(n);
    sortKey(4, numbers, key);
    for (int i = 0; i < n; i++) keyv[i] = key[i];
    std::set<int> uniqueInts2(keyv.begin(), keyv.begin()+4);
    int numUniqueInts2 = uniqueInts2.size();
    REQUIRE(numbers == sortedNumbers);
    REQUIRE(numUniqueInts2 == 4);
    REQUIRE((key[0] == 1 or key[0] == 3));
    REQUIRE((key[1] == 1 or key[1] == 3));
    REQUIRE(key[2] == 0);
    REQUIRE(key[3] == 2);
    delete[] key;

    doubles = {1.1, 2.2, -3.3, -2.2, 2.2, 3.3, -1.1, 0., 2.2, -3.3};
    sortedDoubles = {-3.3, -3.3, -2.2, -1.1, 0., 1.1, 2.2, 2.2, 2.2, 3.3};
    n = doubles.size();
    key = new int[n];
    keyv.resize(n);
    sortKey(n, doubles, key);
    for (int i = 0; i < n; i++) keyv[i] = key[i];
    std::set<int> uniqueDoubles(keyv.begin(), keyv.end());
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
    delete[] key;

    doubles = {1.1, 2.2, -3.3, -2.2, 2.2, 3.3, -1.1, 0., 2.2, -3.3};
    sortedDoubles = {-3.3, -2.2, 1.1, 2.2, 2.2, 3.3, -1.1, 0., 2.2, -3.3};
    n = doubles.size();
    key = new int[n];
    keyv.resize(n);
    sortKey(6, doubles, key);
    for (int i = 0; i < n; i++) keyv[i] = key[i];
    std::set<int> uniqueDoubles2(keyv.begin(), keyv.begin()+6);
    int numUniqueDoubles2 = uniqueDoubles2.size();
    REQUIRE(doubles == sortedDoubles);
    REQUIRE(numUniqueDoubles2 == 6);
    REQUIRE(key[0] == 2);
    REQUIRE(key[1] == 3);
    REQUIRE(key[2] == 0);
    REQUIRE((key[3] == 1 or key[3] == 4));
    REQUIRE((key[4] == 1 or key[4] == 4));
    REQUIRE(key[5] == 5);
    delete[] key;
}
}
