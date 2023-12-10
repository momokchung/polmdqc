#include "catch.hpp"
#include "numeral.h"
#include <string>

namespace polmdqc
{
TEST_CASE("numeral") {
    int number = 42;
    int width = 6;
    std::string paddedNumber;
    paddedNumber = numeral(number, width);
    REQUIRE(paddedNumber == "000042");

    width = 2;
    paddedNumber = numeral(number, width);
    REQUIRE(paddedNumber == "42");

    width = 1;
    paddedNumber = numeral(number, width);
    REQUIRE(paddedNumber == "42");
}
}
