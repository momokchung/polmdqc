#include "catch.hpp"
#include "getnumb.h"
#include <string>

namespace polmdqc
{
TEST_CASE("getnumb") {
    std::string string = "a 42 ";
    int number;
    int next = 0;
    getnumb(string, number, next);
    REQUIRE(number == 0);
    REQUIRE(next == 0);
    next = 1;
    getnumb(string, number, next);
    REQUIRE(number == 42);
    REQUIRE(next == 4);
    string = "  42a";
    next = 2;
    getnumb(string, number, next);
    REQUIRE(number == 42);
    REQUIRE(next == 4);
    string = "    a42   2";
    next = 3;
    getnumb(string, number, next);
    REQUIRE(number == 0);
    REQUIRE(next == 3);
    next = 8;
    getnumb(string, number, next);
    REQUIRE(number == 2);
    REQUIRE(next == 11);
    string = "atom          5    6    O     \"Glycine O\"                    8    15.995    1";
    next = 4;
    getnumb(string, number, next);
    REQUIRE(number == 5);
    REQUIRE(next == 15);
    getnumb(string, number, next);
    REQUIRE(number == 6);
    REQUIRE(next == 20);
    string = "atom          5    6";
    next = 4;
    getnumb(string, number, next);
    REQUIRE(number == 5);
    REQUIRE(next == 15);
    getnumb(string, number, next);
    REQUIRE(number == 6);
    REQUIRE(next == 20);
}
}
