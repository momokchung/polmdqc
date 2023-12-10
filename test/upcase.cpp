#include "catch.hpp"
#include "upcase.h"
#include <string>

namespace polmdqc
{
TEST_CASE("upcase") {
    std::string string = " Hello, World!  ";
    upcase(string);
    REQUIRE(string == " HELLO, WORLD!  ");
}
}
