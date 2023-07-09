#include "catch.hpp"
#include "trimtext.h"
#include <string>

TEST_CASE("trimtext") {
    std::string string = " Hello, World!  \n";
    int trimint = trimtext(string);
    REQUIRE(trimint == 13);
    string = "      ";
    trimint = trimtext(string);
    REQUIRE(trimint == -1);
}
