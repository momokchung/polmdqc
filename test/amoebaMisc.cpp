#include "action.h"
#include "amoebaMisc.h"
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
TEST_CASE("amoebaMisc-1", "[analyze][AMOEBA][water09]") {
    // tests ability to find water09.xyz_2 if input is just water09
    int argc = 3;
    const char* strings[] = {
        "analyze",
        "../../test/testFiles/amoebaMisc/water09",
        "e",
    };
    char** argv = const_cast<char**>(strings);

    analyze(argc, argv);

    REQUIRE(nem == amoebaMisc1::nem);

    COMPARE_REALS(einter, amoebaMisc1::einter, amoebaMisc1::eps);
    COMPARE_REALS(esum, amoebaMisc1::esum, amoebaMisc1::eps);
    COMPARE_REALS(em, amoebaMisc1::em, amoebaMisc1::eps);

    COMPARE_VECTOR(aesum, amoebaMisc1::aesum, amoebaMisc1::eps);
    COMPARE_VECTOR(aem, amoebaMisc1::aem, amoebaMisc1::eps);
}

TEST_CASE("amoebaMisc-2", "[analyze][AMOEBA][alatet_water09]") {
    // tests ability to read concatenated files
    int argc = 3;
    const char* strings[] = {
        "analyze",
        "../../test/testFiles/amoebaMisc/alatet_water09.xyz",
        "e",
    };
    char** argv = const_cast<char**>(strings);

    analyze(argc, argv);

    REQUIRE(nem == amoebaMisc2::nem);

    COMPARE_REALS(einter, amoebaMisc2::einter, amoebaMisc2::eps);
    COMPARE_REALS(esum, amoebaMisc2::esum, amoebaMisc2::eps);
    COMPARE_REALS(em, amoebaMisc2::em, amoebaMisc2::eps);

    COMPARE_VECTOR(aesum, amoebaMisc2::aesum, amoebaMisc2::eps);
    COMPARE_VECTOR(aem, amoebaMisc2::aem, amoebaMisc2::eps);
}

TEST_CASE("amoebaMisc-3", "[testgrad][AMOEBA][alatet_water09]") {
    // tests ability to read concatenated files
    int argc = 5;
    const char* strings[] = {
        "testgrad",
        "../../test/testFiles/amoebaMisc/alatet_water09.xyz",
        "Y",
        "Y",
        "1e-5",
    };
    char** argv = const_cast<char**>(strings);

    testgrad(argc, argv);

    COMPARE_MATRIX(desum, amoebaMisc3::desum, amoebaMisc3::eps1);
    COMPARE_MATRIX(dem, amoebaMisc3::dem, amoebaMisc3::eps1);

    COMPARE_MATRIX(ndesum, amoebaMisc3::desum, amoebaMisc3::eps2);
    COMPARE_MATRIX(ndem, amoebaMisc3::dem, amoebaMisc3::eps2);
}

TEST_CASE("amoebaMisc-5", "[analyze][AMOEBA][water09NSeq]") {
    // tests ability to read non-sequential xyz file
    int argc = 3;
    const char* strings[] = {
        "analyze",
        "../../test/testFiles/amoebaMisc/water09NSeq.xyz",
        "e",
    };
    char** argv = const_cast<char**>(strings);

    analyze(argc, argv);

    REQUIRE(nem == amoebaMisc5::nem);

    COMPARE_REALS(einter, amoebaMisc5::einter, amoebaMisc5::eps);
    COMPARE_REALS(esum, amoebaMisc5::esum, amoebaMisc5::eps);
    COMPARE_REALS(em, amoebaMisc5::em, amoebaMisc5::eps);

    COMPARE_VECTOR(aesum, amoebaMisc5::aesum, amoebaMisc5::eps);
    COMPARE_VECTOR(aem, amoebaMisc5::aem, amoebaMisc5::eps);
}
}
