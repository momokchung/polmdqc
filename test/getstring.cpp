#include "catch.hpp"
#include "getstring.h"
#include <string>

namespace polmdqc
{
TEST_CASE("getstring") {
    std::string string = "atom          1    1    N   lool  \"Glycine lol N\"      \"lol    \"          7    14.003    3";
    std::string text;
    int next = 0;
    getstring(string, text, next);
    REQUIRE(text == "Glycine lol N");
    REQUIRE(next == 49);
    next = 25;
    getstring(string, text, next);
    REQUIRE(text == "Glycine lol N");
    REQUIRE(next == 49);
    next = 30;
    getstring(string, text, next);
    REQUIRE(text == "Glycine lol N");
    REQUIRE(next == 49);
    string = "atom          1    1    N   lool  \"Glycine lol N\"";
    next = 0;
    getstring(string, text, next);
    REQUIRE(text == "Glycine lol N");
    REQUIRE(next == 49);
    next = 25;
    getstring(string, text, next);
    REQUIRE(text == "Glycine lol N");
    REQUIRE(next == 49);
    next = 30;
    getstring(string, text, next);
    REQUIRE(text == "Glycine lol N");
    REQUIRE(next == 49);
}
}
