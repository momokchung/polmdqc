#include "action.h"
#include "analyze.h"
#include "energi.h"
#include "inter.h"
#include "testrt.h"
#include <cmath>

namespace polmdqc
{
TEST_CASE("amoeba-1", "[AMOEBA][water09]") {
    int argc = 3;
    const char* strings[] = {
        "analyze",
        "../../test/testFiles/water09/water09.xyz",
        "e"
    };
    char** argv = const_cast<char**>(strings);

    double einterTest = -5.5415386409028713;
    double emTest = -5.5415386409028713;
    double nemTest = 9;
    double eps = 1e-10;

    analyze(argc, argv);

    COMPARE_REALS(einter, einterTest, eps);
    COMPARE_REALS(em, emTest, eps);
    REQUIRE(nem == nemTest);
}
}
