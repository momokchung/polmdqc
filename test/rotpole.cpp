#include "action.h"
#include "analyz.h"
#include "analyze.h"
#include "deriv.h"
#include "energi.h"
#include "inter.h"
#include "mpole.h"
#include "rotpole.h"
#include "testgrad.h"
#include "testrt.h"

namespace polmdqc
{
TEST_CASE("rotpole-1", "[AMOEBA][axetyp]") {
    // Tests None, Bisector, and Z-then-X
    int argc = 3;
    const char* strings[] = {
        "analyze",
        "../../test/testFiles/rotpole/water09_Na.xyz",
        "e"
    };
    char** argv = const_cast<char**>(strings);

    analyze(argc, argv);

    REQUIRE(nem == rotpole1::nem);

    COMPARE_REALS(einter, rotpole1::einter, rotpole1::eps);
    COMPARE_REALS(esum, rotpole1::esum, rotpole1::eps);
    COMPARE_REALS(em, rotpole1::em, rotpole1::eps);

    COMPARE_VECTOR(aesum, rotpole1::aesum, rotpole1::eps);
    COMPARE_VECTOR(aem, rotpole1::aem, rotpole1::eps);

    COMPARE_MATRIX(rpole, rotpole1::rpole, rotpole1::epsR);
}

TEST_CASE("rotpole-2", "[HIPPO][benzene_ethyne_water]") {
    // Tests Z-Only
    int argc = 3;
    const char* strings[] = {
        "analyze",
        "../../test/testFiles/rotpole/benzene_ethyne_water.xyz",
        "e"
    };
    char** argv = const_cast<char**>(strings);

    analyze(argc, argv);

    REQUIRE(nem == rotpole2::nem);

    COMPARE_REALS(einter, rotpole2::einter, rotpole2::eps);
    COMPARE_REALS(esum, rotpole2::esum, rotpole2::eps);
    COMPARE_REALS(em, rotpole2::em, rotpole2::eps);

    COMPARE_VECTOR(aesum, rotpole2::aesum, rotpole2::eps);
    COMPARE_VECTOR(aem, rotpole2::aem, rotpole2::eps);

    COMPARE_MATRIX(rpole, rotpole2::rpole, rotpole2::epsR);
}

TEST_CASE("rotpole-3", "[HIPPO][ammonia]") {
    // Tests Z-Bisect and 3-Fold
    int argc = 3;
    const char* strings[] = {
        "analyze",
        "../../test/testFiles/rotpole/ammonia.xyz",
        "e"
    };
    char** argv = const_cast<char**>(strings);

    analyze(argc, argv);

    REQUIRE(nem == rotpole3::nem);

    COMPARE_REALS(einter, rotpole3::einter, rotpole3::eps);
    COMPARE_REALS(esum, rotpole3::esum, rotpole3::eps);
    COMPARE_REALS(em, rotpole3::em, rotpole3::eps);

    COMPARE_VECTOR(aesum, rotpole3::aesum, rotpole3::eps);
    COMPARE_VECTOR(aem, rotpole3::aem, rotpole3::eps);

    COMPARE_MATRIX(rpole, rotpole3::rpole, rotpole3::epsR);
}
}
