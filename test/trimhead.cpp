#include "catch.hpp"
#include "trimhead.h"
#include <string>

namespace polmdqc
{
TEST_CASE("trimhead") {
    std::string string = " Hello, World!  \n";
    trimhead(string);
    REQUIRE(string == "Hello, World!  \n");
    string = "      \n";
    trimhead(string);
    REQUIRE(string == "\n");
}
}
