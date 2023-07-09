#include "catch.hpp"
#include "trimhead.h"
#include <string>

TEST_CASE("trimhead") {
    std::string string = " Hello, World!  \n";
    trimhead(string);
    REQUIRE(string == "Hello, World!  \n");
}
