#include "catch.hpp"
#include "torphase.h"

namespace polmdqc
{
TEST_CASE("torphase") {
    int ft[6] = {1,2,0,3,5,4};
    real vt[6] = {1.1,2.2,3.3,4.4,5.5,6.6};
    real st[6] = {280.,-200.,50.,0.,-3000.,3000.};
    torphase(ft,vt,st);
    real expected_vt[6] = {3.3,1.1,2.2,4.4,6.6,5.5};
    real expected_st[6] = {50.,-80.,160.,0.,120.,-120.};
    for (int i = 0; i < 6; i++) {
        REQUIRE(vt[i] == expected_vt[i]);
        REQUIRE(st[i] == expected_st[i]);
    }
}
}
