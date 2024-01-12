#include "action.h"
#include "analyze.h"
#include "energi.h"
#include "inter.h"
#include "testrt.h"
#include <cmath>

namespace polmdqc
{
TEST_CASE("hippo-1", "[HIPPO][water21]") {
    int argc = 3;
    const char* strings[] = {
        "analyze",
        "../../test/testFiles/water21/water21.xyz",
        "e"
    };
    char** argv = const_cast<char**>(strings);

    real einterTest = -7.3522485663420181;
    real emTest = -7.3522485663420181;
    real nemTest = 9;
    real eps = 1e-10;

    analyze(argc, argv);

    COMPARE_REALS(einter, einterTest, eps);
    COMPARE_REALS(em, emTest, eps);
    REQUIRE(nem == nemTest);
}
}
