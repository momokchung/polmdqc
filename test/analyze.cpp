#include "action.h"
#include "analyze.h"
#include "energi.h"
#include "inter.h"
#include "testrt.h"
#include <cmath>

namespace polmdqc
{
TEST_CASE("analyze-1", "[AMOEBA][water09]") {
    int argc = 3;
    const char* strings[] = {
        "analyze",
        "../../test/testFiles/water09/water09.xyz",
        "e"
    };
    char** argv = const_cast<char**>(strings);

    double einterTest = -0.25063530556952029;
    double emTest = -0.25063530556952029;
    double nemTest = 27;
    double eps = 1e-13;

    analyze(argc, argv);

    COMPARE_REALS(einter, einterTest, eps);
    COMPARE_REALS(em, emTest, eps);
    REQUIRE(nem == nemTest);
}

TEST_CASE("analyze-2", "[HIPPO][water21]") {
    int argc = 3;
    const char* strings[] = {
        "analyze",
        "../../test/testFiles/water21/water21.xyz",
        "e"
    };
    char** argv = const_cast<char**>(strings);

    double einterTest = -0.22354977915893826;
    double emTest = -0.22354977915893826;
    double nemTest = 27;
    double eps = 1e-12;

    analyze(argc, argv);

    COMPARE_REALS(einter, einterTest, eps);
    COMPARE_REALS(em, emTest, eps);
    REQUIRE(nem == nemTest);
}
}
