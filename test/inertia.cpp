#include "analyze.h"
#include "atomid.h"
#include "atoms.h"
#include "final.h"
#include "inertia.h"
#include "inertiaT.h"
#include "inform.h"
#include "testrt.h"

namespace polmdqc
{
TEST_CASE("inertia-1", "[inertia][AMOEBA][alatet]") {
    int argc = 3;
    const char* strings[] = {
        "analyze",
        "../../test/testFiles/inertia/alatet.xyz",
        "e",
    };
    char** argv = const_cast<char**>(strings);

    test = true;
    analyze(argc, argv);

    int n = 42;
    inertia(2, n, mass.ptr(), x.ptr(), y.ptr(), z.ptr());

    COMPARE_VECTOR(x, inertia1::x, inertia1::eps);
    COMPARE_VECTOR(y, inertia1::y, inertia1::eps);
    COMPARE_VECTOR(z, inertia1::z, inertia1::eps);

    final();
}
}
