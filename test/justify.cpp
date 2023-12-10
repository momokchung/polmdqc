#include "catch.hpp"
#include "justify.h"
#include <string>

namespace polmdqc
{
TEST_CASE("justify") {
    std::string string = " Hello, World!   ";
    justify(string, 50);
    REQUIRE(string == "                                     Hello, World!");
    string = "  ";
    justify(string, 50);
    REQUIRE(string == "                                                  ");
}
}
