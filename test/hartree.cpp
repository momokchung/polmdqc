#include "hartree.h"
#include "kinetic.h"
#include "nuclear.h"
#include "overlap.h"
#include "testrt.h"
#include "tinkerqm.h"

namespace polmdqc
{
TEST_CASE("tinkerqm-1", "[HartreeFock][water_3-21g]") {
    int argc = 3;
    const char* strings[] = {
        "tinkerqm",
        "../../test/testFiles/water_3-21g/water_3-21g.xyz",
    };
    char** argv = const_cast<char**>(strings);

    tinkerqm(argc, argv);

    COMPARE_MATRIX(overlap::cartS, hartree1::cartS, hartree1::epsS);
    COMPARE_MATRIX(kinetic::cartKE, hartree1::cartKE, hartree1::epsKE);
    COMPARE_MATRIX(nuclear::cartNE, hartree1::cartNE, hartree1::epsNE);
}

TEST_CASE("tinkerqm-2", "[HartreeFock][water_aug-cc-pvtz]") {
    int argc = 3;
    const char* strings[] = {
        "tinkerqm",
        "../../test/testFiles/water_aug-cc-pvtz/water_aug-cc-pvtz.xyz",
    };
    char** argv = const_cast<char**>(strings);

    tinkerqm(argc, argv);

    COMPARE_MATRIX(overlap::sphS, hartree2::sphS, hartree2::epsS);
    COMPARE_MATRIX(kinetic::sphKE, hartree2::sphKE, hartree2::epsKE);
    COMPARE_MATRIX(nuclear::sphNE, hartree2::sphNE, hartree2::epsNE);
}
}
