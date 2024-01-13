#include "action.h"
#include "amoeba.h"
#include "analyz.h"
#include "analyze.h"
#include "deriv.h"
#include "energi.h"
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

    analyze(argc, argv);

    COMPARE_REALS(einter, amoeba1::einter, amoeba1::eps);
    COMPARE_REALS(em, amoeba1::em, amoeba1::eps);
    REQUIRE(nem == amoeba1::nem);
    COMPARE_VECTOR(aem, amoeba1::aem, amoeba1::eps);
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

    testgrad(argc, argv);

    COMPARE_MATRIX(desum, amoeba2::desum, amoeba2::eps1);
    COMPARE_MATRIX(dem, amoeba2::dem, amoeba2::eps1);

    COMPARE_MATRIX(ndesum, amoeba2::desum, amoeba2::eps2);
    COMPARE_MATRIX(ndem, amoeba2::dem, amoeba2::eps2);
}
}
