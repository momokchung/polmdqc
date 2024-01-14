#include "action.h"
#include "analyze.h"
#include "energi.h"
#include "inter.h"
#include "mpole.h"
#include "rotpole.h"
#include "testrt.h"

namespace polmdqc
{
TEST_CASE("rotpole-1", "[AMOEBA][axetyp]") {
    // Tests Bisector, None, Z-Only, and Z-then-X
    int argc = 3;
    const char* strings[] = {
        "analyze",
        "../../test/testFiles/rotpole/axetyp.xyz",
        "e"
    };
    char** argv = const_cast<char**>(strings);

    analyze(argc, argv);

    COMPARE_MATRIX(rpole, rotpole1::rpole, rotpole1::epsR);
    COMPARE_REALS(einter, rotpole1::einter, rotpole1::eps);
    COMPARE_REALS(em, rotpole1::em, rotpole1::eps);
    REQUIRE(nem == rotpole1::nem);
}

TEST_CASE("rotpole-2", "[AMOEBA][lysine_zbisect]") {
    int argc = 3;
    const char* strings[] = {
        "analyze",
        "../../test/testFiles/rotpole/lysine_zbisect.xyz",
        "e"
    };
    char** argv = const_cast<char**>(strings);

    analyze(argc, argv);

    COMPARE_VECTOR(rpole[0], rotpole2::rpole, rotpole2::epsR);
    COMPARE_REALS(em, rotpole2::em, rotpole2::eps);
    REQUIRE(nem == rotpole2::nem);
}

TEST_CASE("rotpole-3", "[AMOEBA][lysine_3fold]") {
    int argc = 3;
    const char* strings[] = {
        "analyze",
        "../../test/testFiles/rotpole/lysine_3fold.xyz",
        "e"
    };
    char** argv = const_cast<char**>(strings);

    analyze(argc, argv);

    COMPARE_VECTOR(rpole[0], rotpole3::rpole, rotpole3::epsR);
    COMPARE_REALS(em, rotpole3::em, rotpole3::eps);
    REQUIRE(nem == rotpole3::nem);
}
}
