#include "action.h"
#include "amoebap.h"
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
TEST_CASE("amoebap-1", "[analyze][AMOEBAPLUS][waterap]") {
    int argc = 3;
    const char* strings[] = {
        "analyze",
        "../../test/testFiles/amoebap/waterap.xyz",
        "e"
    };
    char** argv = const_cast<char**>(strings);

    analyze(argc, argv);

    REQUIRE(nem == amoebap1::nem);

    COMPARE_REALS(einter, amoebap1::einter, amoebap1::eps);
    COMPARE_REALS(esum, amoebap1::esum, amoebap1::eps);
    COMPARE_REALS(em, amoebap1::em, amoebap1::eps);

    COMPARE_VECTOR(aesum, amoebap1::aesum, amoebap1::eps);
    COMPARE_VECTOR(aem, amoebap1::aem, amoebap1::eps);
}

TEST_CASE("amoebap-2", "[testgrad][AMOEBAPLUS][waterap]") {
    int argc = 5;
    const char* strings[] = {
        "testgrad",
        "../../test/testFiles/amoebap/waterap.xyz",
        "Y",
        "Y",
        "1e-5",
    };
    char** argv = const_cast<char**>(strings);

    testgrad(argc, argv);

    COMPARE_MATRIX(desum, amoebap2::desum, amoebap2::eps1);
    COMPARE_MATRIX(dem, amoebap2::dem, amoebap2::eps1);

    COMPARE_MATRIX(ndesum, amoebap2::desum, amoebap2::eps2);
    COMPARE_MATRIX(ndem, amoebap2::dem, amoebap2::eps2);
}
}
