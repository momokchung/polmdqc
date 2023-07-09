#include "catch.hpp"
#include "lowcase.h"
#include <string>

TEST_CASE("lowcase") {
    std::string string = " Hello, World!  ";
    lowcase(string);
    REQUIRE(string == " hello, world!  ");
}
