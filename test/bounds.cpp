#include "action.h"
#include "analyz.h"
#include "analyze.h"
#include "bounds.h"
#include "energi.h"
#include "inter.h"
#include "testrt.h"

namespace polmdqc
{
TEST_CASE("bounds-1", "[analyze][AMOEBA][orthogonal]") {
    // Tests orthogonal
    int argc = 5;
    const char* strings[] = {
        "analyze",
        "-k",
        "../../test/testFiles/bounds/orthogonal.key",
        "../../test/testFiles/bounds/waterbox30.xyz",
        "e"
    };
    char** argv = const_cast<char**>(strings);

    analyze(argc, argv);

    REQUIRE(nem == bounds1::nem);

    COMPARE_REALS(einter, bounds1::einter, bounds1::eps);
    COMPARE_REALS(esum, bounds1::esum, bounds1::eps);
    COMPARE_REALS(em, bounds1::em, bounds1::eps);
}

TEST_CASE("bounds-2", "[analyze][AMOEBA][monoclinic]") {
    // Tests monoclinic
    int argc = 5;
    const char* strings[] = {
        "analyze",
        "-k",
        "../../test/testFiles/bounds/monoclinic.key",
        "../../test/testFiles/bounds/waterbox30.xyz",
        "e"
    };
    char** argv = const_cast<char**>(strings);

    analyze(argc, argv);

    REQUIRE(nem == bounds2::nem);

    COMPARE_REALS(einter, bounds2::einter, bounds2::eps);
    COMPARE_REALS(esum, bounds2::esum, bounds2::eps);
    COMPARE_REALS(em, bounds2::em, bounds2::eps);
}

TEST_CASE("bounds-3", "[analyze][AMOEBA][triclinic]") {
    // Tests triclinic
    int argc = 5;
    const char* strings[] = {
        "analyze",
        "-k",
        "../../test/testFiles/bounds/triclinic.key",
        "../../test/testFiles/bounds/waterbox30.xyz",
        "e"
    };
    char** argv = const_cast<char**>(strings);

    analyze(argc, argv);

    REQUIRE(nem == bounds3::nem);

    COMPARE_REALS(einter, bounds3::einter, bounds3::eps);
    COMPARE_REALS(esum, bounds3::esum, bounds3::eps);
    COMPARE_REALS(em, bounds3::em, bounds3::eps);
}

TEST_CASE("bounds-4", "[analyze][AMOEBA][octahedron]") {
    // Tests octahedron
    int argc = 5;
    const char* strings[] = {
        "analyze",
        "-k",
        "../../test/testFiles/bounds/octahedron.key",
        "../../test/testFiles/bounds/waterbox30.xyz",
        "e"
    };
    char** argv = const_cast<char**>(strings);

    analyze(argc, argv);

    REQUIRE(nem == bounds4::nem);

    COMPARE_REALS(einter, bounds4::einter, bounds4::eps);
    COMPARE_REALS(esum, bounds4::esum, bounds4::eps);
    COMPARE_REALS(em, bounds4::em, bounds4::eps);
}

TEST_CASE("bounds-5", "[analyze][AMOEBA][dodecahedron]") {
    // Tests dodecahedron
    int argc = 5;
    const char* strings[] = {
        "analyze",
        "-k",
        "../../test/testFiles/bounds/dodecahedron.key",
        "../../test/testFiles/bounds/waterbox30.xyz",
        "e"
    };
    char** argv = const_cast<char**>(strings);

    analyze(argc, argv);

    REQUIRE(nem == bounds5::nem);

    COMPARE_REALS(einter, bounds5::einter, bounds5::eps);
    COMPARE_REALS(esum, bounds5::esum, bounds5::eps);
    COMPARE_REALS(em, bounds5::em, bounds5::eps);
}
}
