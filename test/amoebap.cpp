#include "action.h"
#include "analyze.h"
#include "energi.h"
#include "inter.h"
#include "testrt.h"
#include <cmath>

namespace polmdqc
{
TEST_CASE("amoebap-1", "[AMOEBAPLUS][water]") {
    int argc = 3;
    const char* strings[] = {
        "analyze",
        "../../test/testFiles/waterap/waterap.xyz",
        "e"
    };
    char** argv = const_cast<char**>(strings);

    double einterTest = -7.5295245698070445;
    double emTest = -7.5295245698070445;
    double nemTest = 9;
    double eps = 1e-10;

    analyze(argc, argv);

    COMPARE_REALS(einter, einterTest, eps);
    COMPARE_REALS(em, emTest, eps);
    REQUIRE(nem == nemTest);
}
}
