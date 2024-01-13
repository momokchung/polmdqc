#include "action.h"
#include "analyz.h"
#include "analyze.h"
#include "deriv.h"
#include "energi.h"
#include "hippo.h"
#include "inter.h"
#include "testgrad.h"
#include "testrt.h"
#include <cmath>

namespace polmdqc
{
TEST_CASE("hippo-1", "[analyze][HIPPO][water21]") {
    int argc = 3;
    const char* strings[] = {
        "analyze",
        "../../test/testFiles/hippo/water21.xyz",
        "e",
    };
    char** argv = const_cast<char**>(strings);

    analyze(argc, argv);

    COMPARE_REALS(einter, hippo1::einter, hippo1::eps);
    COMPARE_REALS(em, hippo1::em, hippo1::eps);
    REQUIRE(nem == hippo1::nem);
    COMPARE_VECTOR(aem, hippo1::aem, hippo1::eps);
}

TEST_CASE("hippo-2", "[testgrad][HIPPO][water21]") {
    int argc = 5;
    const char* strings[] = {
        "testgrad",
        "../../test/testFiles/hippo/water21.xyz",
        "Y",
        "Y",
        "1e-5",
    };
    char** argv = const_cast<char**>(strings);

    testgrad(argc, argv);

    COMPARE_MATRIX(desum, hippo2::desum, hippo2::eps1);
    COMPARE_MATRIX(dem, hippo2::dem, hippo2::eps1);

    COMPARE_MATRIX(ndesum, hippo2::desum, hippo2::eps2);
    COMPARE_MATRIX(ndem, hippo2::dem, hippo2::eps2);
}
}
