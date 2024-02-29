// Author: Moses KJ Chung
// Year:   2024

#pragma once

namespace polmdqc
{
////////////////////////////////////////////////////////////////////
//                                                                //
//  regular_convex  --  check for local regularity and convexity  //
//                                                                //
////////////////////////////////////////////////////////////////////

void regular_convex(int a, int b, int c, int p, int o, int itest_abcp,
    bool& regular, bool& convex, bool& test_abpo, bool& test_bcpo, bool& test_capo) ;
}
