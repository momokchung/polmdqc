#include "boyscheb.h"
#include "boysref.h"
#include "cheb.h"
#include "inform.h"
#include "initialqm.h"
#include "libfunc.h"
#include "testrt.h"

namespace polmdqc
{
TEST_CASE("boyscheb-1", "[Chebyshev]") {
    int argc = 0;
    const char* strings[] = {};
    char** argv = const_cast<char**>(strings);
    test = true;

    // allocate dynamic arrays
    double* FmCheb;
    double* FmRef;
    FmCheb = new double[chebmmax+1];
    FmRef = new double[chebmmax+1];

    // set number of intervals
    const int interval = 840;

    // set shift and epsilon
    double shft = 0.02857142857; // shift = 0.2/7
    double abseps = 1e-15;
    double releps = 1e-13;

    // initialize
    initialqm(argc, argv);
    double maxabsdiff = 0.;
    double maxreldiff = 0.;

    for (int i = 0; i < interval; i++) {
        double midpoint = i * chebdlta + chebdlta2;
        double t = midpoint + shft;
        boyscheb(FmCheb,chebmmax,t);
        boysref(FmRef,chebmmax,t);
        for (int m = 0; m <= chebmmax; m++) {
            double fc = FmCheb[m];
            double fr = FmRef[m];
            double absdiff = REAL_ABS(fc - fr);
            double reldiff = absdiff / fr;
            // printf("i m %d %d\n", i, m);
            // printf("absdiff reldiff %25.16e%25.16e\n", absdiff,reldiff);
            REQUIRE(absdiff < abseps);
            REQUIRE(reldiff < releps);
            if (absdiff > maxabsdiff) maxabsdiff = absdiff;
            if (reldiff > maxreldiff) maxreldiff = reldiff;
        }
    }

    // printf("maxabsdiff, maxreldiff: %25.16e%25.16e\n", maxabsdiff, maxreldiff);

    // deallocate dynamic arrays
    delete[] FmCheb;
    delete[] FmRef;
}
}
