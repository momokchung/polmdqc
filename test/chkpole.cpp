#include "action.h"
#include "analyz.h"
#include "analyze.h"
#include "chkpole.h"
#include "deriv.h"
#include "energi.h"
#include "final.h"
#include "inform.h"
#include "inter.h"
#include "testgrad.h"
#include "testrt.h"
#include <cmath>

namespace polmdqc
{
TEST_CASE("chkpole-1", "[analyze][AMOEBA][alatet]") {
    int argc = 3;
    const char* strings[] = {
        "analyze",
        "../../test/testFiles/chkpole/alatet.xyz",
        "e",
    };
    char** argv = const_cast<char**>(strings);

    test = true;
    analyze(argc, argv);

    REQUIRE(nem == chkpole1::nem);

    COMPARE_REALS(einter, chkpole1::einter, chkpole1::eps);
    COMPARE_REALS(esum, chkpole1::esum, chkpole1::eps);
    COMPARE_REALS(em, chkpole1::em, chkpole1::eps);

    COMPARE_VECTOR(aesum, chkpole1::aesum, chkpole1::eps);
    COMPARE_VECTOR(aem, chkpole1::aem, chkpole1::eps);

    final();
}

TEST_CASE("chkpole-2", "[testgrad][AMOEBA][alatet]") {
    int argc = 5;
    const char* strings[] = {
        "testgrad",
        "../../test/testFiles/chkpole/alatet.xyz",
        "Y",
        "Y",
        "1e-5",
    };
    char** argv = const_cast<char**>(strings);

    test = true;
    testgrad(argc, argv);

    COMPARE_VECTOR(desum, chkpole2::desum, chkpole2::eps1);
    COMPARE_VECTOR(dem, chkpole2::dem, chkpole2::eps1);

    COMPARE_VECTOR(ndesum, chkpole2::desum, chkpole2::eps2);
    COMPARE_VECTOR(ndem, chkpole2::dem, chkpole2::eps2);

    final();
}
}
