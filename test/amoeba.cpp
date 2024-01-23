#include "action.h"
#include "amoeba.h"
#include "analyz.h"
#include "analyze.h"
#include "deriv.h"
#include "energi.h"
#include "final.h"
#include "inter.h"
#include "testgrad.h"
#include "testrt.h"
#include <cmath>

namespace polmdqc
{
TEST_CASE("amoeba-1", "[analyze][AMOEBA][water09]") {
    int argc = 3;
    const char* strings[] = {
        "analyze",
        "../../test/testFiles/amoeba/water09.xyz",
        "e",
    };
    char** argv = const_cast<char**>(strings);

    analyze(argc, argv, true);

    REQUIRE(nem == amoeba1::nem);

    COMPARE_REALS(einter, amoeba1::einter, amoeba1::eps);
    COMPARE_REALS(esum, amoeba1::esum, amoeba1::eps);
    COMPARE_REALS(em, amoeba1::em, amoeba1::eps);

    COMPARE_VECTOR(aesum, amoeba1::aesum, amoeba1::eps);
    COMPARE_VECTOR(aem, amoeba1::aem, amoeba1::eps);

    final();
}

TEST_CASE("amoeba-2", "[testgrad][AMOEBA][water09]") {
    int argc = 5;
    const char* strings[] = {
        "testgrad",
        "../../test/testFiles/amoeba/water09.xyz",
        "Y",
        "Y",
        "1e-5",
    };
    char** argv = const_cast<char**>(strings);

    testgrad(argc, argv, true);

    COMPARE_VECTOR(desum, amoeba2::desum, amoeba2::eps1);
    COMPARE_VECTOR(dem, amoeba2::dem, amoeba2::eps1);

    COMPARE_VECTOR(ndesum, amoeba2::desum, amoeba2::eps2);
    COMPARE_VECTOR(ndem, amoeba2::dem, amoeba2::eps2);

    final();
}
}
