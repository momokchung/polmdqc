#include "action.h"
#include "active.h"
#include "analyz.h"
#include "analyze.h"
#include "atoms.h"
#include "deriv.h"
#include "energi.h"
#include "final.h"
#include "inform.h"
#include "inter.h"
#include "testgrad.h"
#include "testrt.h"
#include "usage.h"
#include <cmath>

namespace polmdqc
{
TEST_CASE("active-1", "[analyze][AMOEBA][water09_Na_Cls]") {
    // tests active keyword
    int argc = 5;
    const char* strings[] = {
        "analyze",
        "-k",
        "../../test/testFiles/active/active.key",
        "../../test/testFiles/active/water09_Na_Cls.xyz",
        "e",
    };
    char** argv = const_cast<char**>(strings);

    test = true;
    analyze(argc, argv);

    REQUIRE(nem == active1::nem);

    COMPARE_REALS(einter, active1::einter, active1::eps);
    COMPARE_REALS(esum, active1::esum, active1::eps);
    COMPARE_REALS(em, active1::em, active1::eps);

    COMPARE_VECTOR(aesum, active1::aesum, active1::eps);
    COMPARE_VECTOR(aem, active1::aem, active1::eps);

    REQUIRE(nuse == active1::nuse);
    for (int i = 0; i < n; i++) {
        REQUIRE(iuse[i] == active1::iuse[i]);
        REQUIRE(use[i+1] == active1::use[i+1]);
    }

    final();
}

TEST_CASE("active-2", "[analyze][AMOEBA][water09_Na_Cls]") {
    // tests inactive keyword
    int argc = 5;
    const char* strings[] = {
        "analyze",
        "-k",
        "../../test/testFiles/active/inactive.key",
        "../../test/testFiles/active/water09_Na_Cls.xyz",
        "e",
    };
    char** argv = const_cast<char**>(strings);

    test = true;
    analyze(argc, argv);

    REQUIRE(nem == active2::nem);

    COMPARE_REALS(einter, active2::einter, active2::eps);
    COMPARE_REALS(esum, active2::esum, active2::eps);
    COMPARE_REALS(em, active2::em, active2::eps);

    COMPARE_VECTOR(aesum, active2::aesum, active2::eps);
    COMPARE_VECTOR(aem, active2::aem, active2::eps);

    REQUIRE(nuse == active2::nuse);
    for (int i = 0; i < n; i++) {
        REQUIRE(iuse[i] == active2::iuse[i]);
        REQUIRE(use[i+1] == active2::use[i+1]);
    }

    final();
}

TEST_CASE("active-3", "[analyze][AMOEBA][waterbox30]") {
    // tests active-sphere keyword
    int argc = 3;
    const char* strings[] = {
        "analyze",
        "../../test/testFiles/active/waterbox30.xyz",
        "e",
    };
    char** argv = const_cast<char**>(strings);

    test = true;
    analyze(argc, argv);

    REQUIRE(nem == active3::nem);

    COMPARE_REALS(einter, active3::einter, active3::eps);
    COMPARE_REALS(esum, active3::esum, active3::eps);
    COMPARE_REALS(em, active3::em, active3::eps);

    REQUIRE(nuse == active3::nuse);

    final();
}
}
