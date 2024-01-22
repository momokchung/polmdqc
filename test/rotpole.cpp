#include "action.h"
#include "analyz.h"
#include "analyze.h"
#include "deriv.h"
#include "energi.h"
#include "final.h"
#include "inter.h"
#include "mpole.h"
#include "rotpole.h"
#include "testgrad.h"
#include "testrt.h"

namespace polmdqc
{
TEST_CASE("rotpole-1", "[analyze][AMOEBA][axetyp]") {
    // Tests None, Bisector, and Z-then-X
    int argc = 3;
    const char* strings[] = {
        "analyze",
        "../../test/testFiles/rotpole/water09_Na.xyz",
        "e"
    };
    char** argv = const_cast<char**>(strings);

    analyze(argc, argv, true);

    REQUIRE(nem == rotpole1::nem);

    COMPARE_REALS(einter, rotpole1::einter, rotpole1::eps);
    COMPARE_REALS(esum, rotpole1::esum, rotpole1::eps);
    COMPARE_REALS(em, rotpole1::em, rotpole1::eps);

    COMPARE_VECTOR(aesum, rotpole1::aesum, rotpole1::eps);
    COMPARE_VECTOR(aem, rotpole1::aem, rotpole1::eps);

    COMPARE_ARRAY2D(rpole, rotpole1::rpole, rotpole1::epsR);

    final();
}

TEST_CASE("rotpole-2", "[analyze][HIPPO][benzene_ethyne_water]") {
    // Tests Z-Only
    int argc = 3;
    const char* strings[] = {
        "analyze",
        "../../test/testFiles/rotpole/benzene_ethyne_water.xyz",
        "e"
    };
    char** argv = const_cast<char**>(strings);

    analyze(argc, argv, true);

    REQUIRE(nem == rotpole2::nem);

    COMPARE_REALS(einter, rotpole2::einter, rotpole2::eps);
    COMPARE_REALS(esum, rotpole2::esum, rotpole2::eps);
    COMPARE_REALS(em, rotpole2::em, rotpole2::eps);

    COMPARE_VECTOR(aesum, rotpole2::aesum, rotpole2::eps);
    COMPARE_VECTOR(aem, rotpole2::aem, rotpole2::eps);

    COMPARE_ARRAY2D(rpole, rotpole2::rpole, rotpole2::epsR);

    final();
}

TEST_CASE("rotpole-3", "[analyze][HIPPO][ammonia]") {
    // Tests Z-Bisect and 3-Fold
    int argc = 3;
    const char* strings[] = {
        "analyze",
        "../../test/testFiles/rotpole/ammonia.xyz",
        "e"
    };
    char** argv = const_cast<char**>(strings);

    analyze(argc, argv, true);

    REQUIRE(nem == rotpole3::nem);

    COMPARE_REALS(einter, rotpole3::einter, rotpole3::eps);
    COMPARE_REALS(esum, rotpole3::esum, rotpole3::eps);
    COMPARE_REALS(em, rotpole3::em, rotpole3::eps);

    COMPARE_VECTOR(aesum, rotpole3::aesum, rotpole3::eps);
    COMPARE_VECTOR(aem, rotpole3::aem, rotpole3::eps);

    COMPARE_ARRAY2D(rpole, rotpole3::rpole, rotpole3::epsR);

    final();
}

TEST_CASE("rotpole-4", "[testgrad][AMOEBA][axetyp]") {
    // Tests None, Bisector, and Z-then-X
    int argc = 5;
    const char* strings[] = {
        "testgrad",
        "../../test/testFiles/rotpole/water09_Na.xyz",
        "Y",
        "Y",
        "1e-5",
    };
    char** argv = const_cast<char**>(strings);

    testgrad(argc, argv, true);

    COMPARE_ARRAY2D(desum, rotpole4::desum, rotpole4::eps1);
    COMPARE_ARRAY2D(dem, rotpole4::dem, rotpole4::eps1);

    COMPARE_ARRAY2D(ndesum, rotpole4::desum, rotpole4::eps2);
    COMPARE_ARRAY2D(ndem, rotpole4::dem, rotpole4::eps2);

    final();
}

TEST_CASE("rotpole-5", "[testgrad][HIPPO][benzene_ethyne_water]") {
    // Tests Z-Only
    int argc = 5;
    const char* strings[] = {
        "testgrad",
        "../../test/testFiles/rotpole/benzene_ethyne_water.xyz",
        "Y",
        "Y",
        "1e-5",
    };
    char** argv = const_cast<char**>(strings);

    testgrad(argc, argv, true);

    COMPARE_ARRAY2D(desum, rotpole5::desum, rotpole5::eps1);
    COMPARE_ARRAY2D(dem, rotpole5::dem, rotpole5::eps1);

    COMPARE_ARRAY2D(ndesum, rotpole5::desum, rotpole5::eps2);
    COMPARE_ARRAY2D(ndem, rotpole5::dem, rotpole5::eps2);

    final();
}

TEST_CASE("rotpole-6", "[analyze][HIPPO][ammonia]") {
    // Tests Z-Bisect and 3-Fold
    int argc = 5;
    const char* strings[] = {
        "testgrad",
        "../../test/testFiles/rotpole/ammonia.xyz",
        "Y",
        "Y",
        "1e-5",
    };
    char** argv = const_cast<char**>(strings);

    testgrad(argc, argv, true);

    COMPARE_ARRAY2D(desum, rotpole6::desum, rotpole6::eps1);
    COMPARE_ARRAY2D(dem, rotpole6::dem, rotpole6::eps1);

    COMPARE_ARRAY2D(ndesum, rotpole6::desum, rotpole6::eps2);
    COMPARE_ARRAY2D(ndem, rotpole6::dem, rotpole6::eps2);

    final();
}
}
