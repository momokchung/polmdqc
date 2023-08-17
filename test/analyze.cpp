#include "action.h"
#include "analyze.h"
#include "catch.hpp"
#include "energi.h"
#include "inter.h"
#include <cmath>

TEST_CASE("analyze-1", "[AMOEBA][water09]") {
    int argc = 3;
    const char* strings[] = {
        "analyze",
        "../../test/testFiles/water09/water09.xyz",
        "e"
    };
    char** argv = const_cast<char**>(strings);

    double einterTest = -0.25063531;
    double emTest = -0.25063531;
    double nemTest = 27;

    analyze(argc, argv);

    double einterDiff = einterTest - einter;
    double emDiff = emTest - em;
    int nemDiff = nemTest - nem;

    double eps = 1e-8;
    REQUIRE(std::abs(einterDiff) < eps);
    REQUIRE(std::abs(emDiff) < eps);
    REQUIRE(nemDiff == 0);
}

TEST_CASE("analyze-2", "[HIPPO][water21]") {
    int argc = 3;
    const char* strings[] = {
        "analyze",
        "../../test/testFiles/water21/water21.xyz",
        "e"
    };
    char** argv = const_cast<char**>(strings);

    double einterTest = -0.22354978;
    double emTest = -0.22354978;
    double nemTest = 27;

    analyze(argc, argv);

    double einterDiff = einterTest - einter;
    double emDiff = emTest - em;
    int nemDiff = nemTest - nem;

    double eps = 1e-8;
    REQUIRE(std::abs(einterDiff) < eps);
    REQUIRE(std::abs(emDiff) < eps);
    REQUIRE(nemDiff == 0);
}
