#include "action.h"
#include "analyz.h"
#include "analyze.h"
#include "atoms.h"
#include "deriv.h"
#include "energi.h"
#include "final.h"
#include "hippo.h"
#include "inform.h"
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

    test = true;
    analyze(argc, argv);

    REQUIRE(nem == hippo1::nem);

    COMPARE_REALS(einter, hippo1::einter, hippo1::eps);
    COMPARE_REALS(esum, hippo1::esum, hippo1::eps);
    COMPARE_REALS(em, hippo1::em, hippo1::eps);

    COMPARE_VECTOR(aesum, hippo1::aesum, hippo1::eps);
    COMPARE_VECTOR(aem, hippo1::aem, hippo1::eps);

    final();
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

    test = true;
    testgrad(argc, argv);

    COMPARE_VECTOR(desum, hippo2::desum, hippo2::eps1);
    COMPARE_VECTOR(dem, hippo2::dem, hippo2::eps1);

    COMPARE_VECTOR(ndesum, hippo2::desum, hippo2::eps2);
    COMPARE_VECTOR(ndem, hippo2::dem, hippo2::eps2);

    final();
}

TEST_CASE("hippo-3", "[testgrad][HIPPO][water21]") {
    int argc = 4;
    const char* strings[] = {
        "testgrad",
        "../../test/testFiles/hippo/water21.xyz",
        "Y",
        "N",
    };
    char** argv = const_cast<char**>(strings);

    test = true;
    testgrad(argc, argv);

    COMPARE_REALS(esum, hippo1::esum, hippo1::eps);
    COMPARE_REALS(em, hippo1::em, hippo1::eps);

    COMPARE_VECTOR(desum, hippo2::desum, hippo2::eps1);
    COMPARE_VECTOR(dem, hippo2::dem, hippo2::eps1);

    final();
}
}
