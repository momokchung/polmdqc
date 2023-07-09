#include "catch.hpp"
#include "gettext.h"
#include <string>

TEST_CASE("gettext") {
    std::string string = " Example  string!  ";
    std::string text;
    int next = 0;
    gettext(string, text, next);
    REQUIRE(text == "Example");
    REQUIRE(next == 8);
    gettext(string, text, next);
    REQUIRE(text == "string!");
    REQUIRE(next == 17);
    string = "   \t \n";
    next = 0;
    gettext(string, text, next);
    REQUIRE(text == "");
    REQUIRE(next == 5);
}
