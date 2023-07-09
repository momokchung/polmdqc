#include "catch.hpp"
#include "upcase.h"
#include <string>

TEST_CASE("upcase") {
    std::string string = " Hello, World!  ";
    upcase(string);
    REQUIRE(string == " HELLO, WORLD!  ");
}
