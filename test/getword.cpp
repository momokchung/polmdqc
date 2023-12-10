#include "catch.hpp"
#include "getword.h"
#include <string>

namespace polmdqc
{
TEST_CASE("getword") {
    std::string string = "  Hel,,lo3 4There!!  \n";
    std::string word;
    int next = 0;
    getword(string, word, next);
    REQUIRE(word == "Hel");
    REQUIRE(next == 5);
    getword(string, word, next);
    REQUIRE(word == "lo3");
    REQUIRE(next == 10);
    getword(string, word, next);
    REQUIRE(word == "There!!");
    REQUIRE(next == 19);
    string = "    ";
    next = 0;
    getword(string, word, next);
    REQUIRE(word == "");
    REQUIRE(next == 0);
    string = "Hello";
    next = 0;
    getword(string, word, next);
    REQUIRE(word == "Hello");
    REQUIRE(next == 5);
}
}
