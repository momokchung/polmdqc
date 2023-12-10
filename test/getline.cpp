#include "catch.hpp"
#include "getline.h"
#include <string>

namespace polmdqc
{
TEST_CASE("getline") {
    std::string string = " Hello, World!  \n";
    getline(string);
    REQUIRE(string == "Hello, World!");
    string = "   \t \n";
    getline(string);
    REQUIRE(string == "");
}
}
