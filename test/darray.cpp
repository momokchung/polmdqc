#include "catch.hpp"
#include "darray.h"

namespace polmdqc
{
TEST_CASE("darray-1", "[MDQCArray]") {
    MDQCArray<int> mdqcv1;
    REQUIRE(mdqcv1.size() == 0);

    mdqcv1.allocate(10);
    REQUIRE(mdqcv1.size() == 10);

    mdqcv1.allocate(5);
    REQUIRE(mdqcv1.size() == 5);

    mdqcv1[0] = 0;
    mdqcv1[1] = 1;
    mdqcv1[2] = 2;
    mdqcv1[3] = 3;
    mdqcv1[4] = 4;
    int mdqcv1_0 = mdqcv1[0];
    int mdqcv1_1 = mdqcv1[1];
    int mdqcv1_2 = mdqcv1[2];
    int mdqcv1_3 = mdqcv1[3];
    int mdqcv1_4 = mdqcv1[4];
    REQUIRE(mdqcv1[0] == 0);
    REQUIRE(mdqcv1[1] == 1);
    REQUIRE(mdqcv1[2] == 2);
    REQUIRE(mdqcv1[3] == 3);
    REQUIRE(mdqcv1[4] == 4);
    REQUIRE(mdqcv1_0 == 0);
    REQUIRE(mdqcv1_1 == 1);
    REQUIRE(mdqcv1_2 == 2);
    REQUIRE(mdqcv1_3 == 3);
    REQUIRE(mdqcv1_4 == 4);

    int* ptr = mdqcv1.ptr();
    ptr[0] = 5;
    ptr[1] = 6;
    ptr[2] = 7;
    ptr[3] = 8;
    ptr[4] = 9;
    mdqcv1_0 = ptr[0];
    mdqcv1_1 = ptr[1];
    mdqcv1_2 = ptr[2];
    mdqcv1_3 = ptr[3];
    mdqcv1_4 = ptr[4];
    REQUIRE(mdqcv1[0] == 5);
    REQUIRE(mdqcv1[1] == 6);
    REQUIRE(mdqcv1[2] == 7);
    REQUIRE(mdqcv1[3] == 8);
    REQUIRE(mdqcv1[4] == 9);
    REQUIRE(mdqcv1_0 == 5);
    REQUIRE(mdqcv1_1 == 6);
    REQUIRE(mdqcv1_2 == 7);
    REQUIRE(mdqcv1_3 == 8);
    REQUIRE(mdqcv1_4 == 9);

    mdqcv1[0] = 0;
    mdqcv1[1] = 1;
    mdqcv1[2] = 2;
    mdqcv1[3] = 3;
    mdqcv1[4] = 4;
    ptr = &mdqcv1[2];
    ptr[0] = 5;
    ptr[1] = 6;
    mdqcv1_0 = ptr[0];
    mdqcv1_1 = ptr[1];
    REQUIRE(mdqcv1[0] == 0);
    REQUIRE(mdqcv1[1] == 1);
    REQUIRE(mdqcv1[2] == 5);
    REQUIRE(mdqcv1[3] == 6);
    REQUIRE(mdqcv1[4] == 4);
    REQUIRE(mdqcv1_0 == 5);
    REQUIRE(mdqcv1_1 == 6);

    mdqcv1.deallocate();
    REQUIRE(mdqcv1.size() == 0);
}

TEST_CASE("darray-2", "[MDQCArray2D]") {
    MDQCArray2D<int,3> mdqcv2;
    REQUIRE(mdqcv2.size() == 0);
    REQUIRE(mdqcv2.size2() == 0);

    mdqcv2.allocate(5);
    REQUIRE(mdqcv2.size() == 5);
    REQUIRE(mdqcv2.size2() == 3);

    mdqcv2.allocate(2);
    REQUIRE(mdqcv2.size() == 2);
    REQUIRE(mdqcv2.size2() == 3);

    mdqcv2[0][0] = 0;
    mdqcv2[0][1] = 1;
    mdqcv2[0][2] = 2;
    mdqcv2[1][0] = 3;
    mdqcv2[1][1] = 4;
    mdqcv2[1][2] = 5;
    int mdqcv2_00 = mdqcv2[0][0];
    int mdqcv2_01 = mdqcv2[0][1];
    int mdqcv2_02 = mdqcv2[0][2];
    int mdqcv2_10 = mdqcv2[1][0];
    int mdqcv2_11 = mdqcv2[1][1];
    int mdqcv2_12 = mdqcv2[1][2];
    REQUIRE(mdqcv2[0][0] == 0);
    REQUIRE(mdqcv2[0][1] == 1);
    REQUIRE(mdqcv2[0][2] == 2);
    REQUIRE(mdqcv2[1][0] == 3);
    REQUIRE(mdqcv2[1][1] == 4);
    REQUIRE(mdqcv2[1][2] == 5);
    REQUIRE(mdqcv2_00 == 0);
    REQUIRE(mdqcv2_01 == 1);
    REQUIRE(mdqcv2_02 == 2);
    REQUIRE(mdqcv2_10 == 3);
    REQUIRE(mdqcv2_11 == 4);
    REQUIRE(mdqcv2_12 == 5);

    int (*ptr)[3] = mdqcv2.ptr();
    ptr[0][0] = 4;
    ptr[0][1] = 5;
    ptr[0][2] = 6;
    ptr[1][0] = 7;
    ptr[1][1] = 8;
    ptr[1][2] = 9;
    mdqcv2_00 = ptr[0][0];
    mdqcv2_01 = ptr[0][1];
    mdqcv2_02 = ptr[0][2];
    mdqcv2_10 = ptr[1][0];
    mdqcv2_11 = ptr[1][1];
    mdqcv2_12 = ptr[1][2];
    REQUIRE(mdqcv2[0][0] == 4);
    REQUIRE(mdqcv2[0][1] == 5);
    REQUIRE(mdqcv2[0][2] == 6);
    REQUIRE(mdqcv2[1][0] == 7);
    REQUIRE(mdqcv2[1][1] == 8);
    REQUIRE(mdqcv2[1][2] == 9);
    REQUIRE(mdqcv2_00 == 4);
    REQUIRE(mdqcv2_01 == 5);
    REQUIRE(mdqcv2_02 == 6);
    REQUIRE(mdqcv2_10 == 7);
    REQUIRE(mdqcv2_11 == 8);
    REQUIRE(mdqcv2_12 == 9);

    mdqcv2.deallocate();
    REQUIRE(mdqcv2.size() == 0);
    REQUIRE(mdqcv2.size2() == 0);
}

TEST_CASE("darray-3", "[MDQCArray3D]") {
    MDQCArray3D<int,2,3> mdqcv3;
    REQUIRE(mdqcv3.size() == 0);
    REQUIRE(mdqcv3.size2() == 0);
    REQUIRE(mdqcv3.size3() == 0);

    mdqcv3.allocate(5);
    REQUIRE(mdqcv3.size() == 5);
    REQUIRE(mdqcv3.size2() == 2);
    REQUIRE(mdqcv3.size3() == 3);

    mdqcv3.allocate(1);
    REQUIRE(mdqcv3.size() == 1);
    REQUIRE(mdqcv3.size2() == 2);
    REQUIRE(mdqcv3.size3() == 3);

    mdqcv3[0][0][0] = 0;
    mdqcv3[0][0][1] = 1;
    mdqcv3[0][0][2] = 2;
    mdqcv3[0][1][0] = 3;
    mdqcv3[0][1][1] = 4;
    mdqcv3[0][1][2] = 5;
    int mdqcv3_000 = mdqcv3[0][0][0];
    int mdqcv3_001 = mdqcv3[0][0][1];
    int mdqcv3_002 = mdqcv3[0][0][2];
    int mdqcv3_010 = mdqcv3[0][1][0];
    int mdqcv3_011 = mdqcv3[0][1][1];
    int mdqcv3_012 = mdqcv3[0][1][2];
    REQUIRE(mdqcv3[0][0][0] == 0);
    REQUIRE(mdqcv3[0][0][1] == 1);
    REQUIRE(mdqcv3[0][0][2] == 2);
    REQUIRE(mdqcv3[0][1][0] == 3);
    REQUIRE(mdqcv3[0][1][1] == 4);
    REQUIRE(mdqcv3[0][1][2] == 5);
    REQUIRE(mdqcv3_000 == 0);
    REQUIRE(mdqcv3_001 == 1);
    REQUIRE(mdqcv3_002 == 2);
    REQUIRE(mdqcv3_010 == 3);
    REQUIRE(mdqcv3_011 == 4);
    REQUIRE(mdqcv3_012 == 5);

    int (*ptr)[2][3] = mdqcv3.ptr();
    ptr[0][0][0] = 4;
    ptr[0][0][1] = 5;
    ptr[0][0][2] = 6;
    ptr[0][1][0] = 7;
    ptr[0][1][1] = 8;
    ptr[0][1][2] = 9;
    mdqcv3_000 = ptr[0][0][0];
    mdqcv3_001 = ptr[0][0][1];
    mdqcv3_002 = ptr[0][0][2];
    mdqcv3_010 = ptr[0][1][0];
    mdqcv3_011 = ptr[0][1][1];
    mdqcv3_012 = ptr[0][1][2];
    REQUIRE(mdqcv3[0][0][0] == 4);
    REQUIRE(mdqcv3[0][0][1] == 5);
    REQUIRE(mdqcv3[0][0][2] == 6);
    REQUIRE(mdqcv3[0][1][0] == 7);
    REQUIRE(mdqcv3[0][1][1] == 8);
    REQUIRE(mdqcv3[0][1][2] == 9);
    REQUIRE(mdqcv3_000 == 4);
    REQUIRE(mdqcv3_001 == 5);
    REQUIRE(mdqcv3_002 == 6);
    REQUIRE(mdqcv3_010 == 7);
    REQUIRE(mdqcv3_011 == 8);
    REQUIRE(mdqcv3_012 == 9);

    mdqcv3.deallocate();
    REQUIRE(mdqcv3.size() == 0);
    REQUIRE(mdqcv3.size2() == 0);
    REQUIRE(mdqcv3.size3() == 0);
}

TEST_CASE("darray-4", "[MDQCVector2D]") {
    MDQCVector2D<int> mdqcv2;
    REQUIRE(mdqcv2.size() == 0);
    REQUIRE(mdqcv2.size2() == 0);

    mdqcv2.allocate(5,3);
    REQUIRE(mdqcv2.size() == 5);
    REQUIRE(mdqcv2.size2() == 3);

    mdqcv2.allocate(2,3);
    REQUIRE(mdqcv2.size() == 2);
    REQUIRE(mdqcv2.size2() == 3);

    mdqcv2[0][0] = 0;
    mdqcv2[0][1] = 1;
    mdqcv2[0][2] = 2;
    mdqcv2[1][0] = 3;
    mdqcv2[1][1] = 4;
    mdqcv2[1][2] = 5;
    int mdqcv2_00 = mdqcv2[0][0];
    int mdqcv2_01 = mdqcv2[0][1];
    int mdqcv2_02 = mdqcv2[0][2];
    int mdqcv2_10 = mdqcv2[1][0];
    int mdqcv2_11 = mdqcv2[1][1];
    int mdqcv2_12 = mdqcv2[1][2];
    REQUIRE(mdqcv2[0][0] == 0);
    REQUIRE(mdqcv2[0][1] == 1);
    REQUIRE(mdqcv2[0][2] == 2);
    REQUIRE(mdqcv2[1][0] == 3);
    REQUIRE(mdqcv2[1][1] == 4);
    REQUIRE(mdqcv2[1][2] == 5);
    REQUIRE(mdqcv2_00 == 0);
    REQUIRE(mdqcv2_01 == 1);
    REQUIRE(mdqcv2_02 == 2);
    REQUIRE(mdqcv2_10 == 3);
    REQUIRE(mdqcv2_11 == 4);
    REQUIRE(mdqcv2_12 == 5);

    int** ptr = mdqcv2.ptr();
    ptr[0][0] = 4;
    ptr[0][1] = 5;
    ptr[0][2] = 6;
    ptr[1][0] = 7;
    ptr[1][1] = 8;
    ptr[1][2] = 9;
    mdqcv2_00 = ptr[0][0];
    mdqcv2_01 = ptr[0][1];
    mdqcv2_02 = ptr[0][2];
    mdqcv2_10 = ptr[1][0];
    mdqcv2_11 = ptr[1][1];
    mdqcv2_12 = ptr[1][2];
    REQUIRE(mdqcv2[0][0] == 4);
    REQUIRE(mdqcv2[0][1] == 5);
    REQUIRE(mdqcv2[0][2] == 6);
    REQUIRE(mdqcv2[1][0] == 7);
    REQUIRE(mdqcv2[1][1] == 8);
    REQUIRE(mdqcv2[1][2] == 9);
    REQUIRE(mdqcv2_00 == 4);
    REQUIRE(mdqcv2_01 == 5);
    REQUIRE(mdqcv2_02 == 6);
    REQUIRE(mdqcv2_10 == 7);
    REQUIRE(mdqcv2_11 == 8);
    REQUIRE(mdqcv2_12 == 9);

    mdqcv2.deallocate();
    REQUIRE(mdqcv2.size() == 0);
    REQUIRE(mdqcv2.size2() == 0);
}

TEST_CASE("darray-5", "[MDQCVector3D]") {
    MDQCVector3D<int> mdqcv3;
    REQUIRE(mdqcv3.size() == 0);
    REQUIRE(mdqcv3.size2() == 0);
    REQUIRE(mdqcv3.size3() == 0);

    mdqcv3.allocate(5,3,4);
    REQUIRE(mdqcv3.size() == 5);
    REQUIRE(mdqcv3.size2() == 3);
    REQUIRE(mdqcv3.size3() == 4);

    mdqcv3.allocate(1,2,3);
    REQUIRE(mdqcv3.size() == 1);
    REQUIRE(mdqcv3.size2() == 2);
    REQUIRE(mdqcv3.size3() == 3);

    mdqcv3[0][0][0] = 0;
    mdqcv3[0][0][1] = 1;
    mdqcv3[0][0][2] = 2;
    mdqcv3[0][1][0] = 3;
    mdqcv3[0][1][1] = 4;
    mdqcv3[0][1][2] = 5;
    int mdqcv3_000 = mdqcv3[0][0][0];
    int mdqcv3_001 = mdqcv3[0][0][1];
    int mdqcv3_002 = mdqcv3[0][0][2];
    int mdqcv3_010 = mdqcv3[0][1][0];
    int mdqcv3_011 = mdqcv3[0][1][1];
    int mdqcv3_012 = mdqcv3[0][1][2];
    REQUIRE(mdqcv3[0][0][0] == 0);
    REQUIRE(mdqcv3[0][0][1] == 1);
    REQUIRE(mdqcv3[0][0][2] == 2);
    REQUIRE(mdqcv3[0][1][0] == 3);
    REQUIRE(mdqcv3[0][1][1] == 4);
    REQUIRE(mdqcv3[0][1][2] == 5);
    REQUIRE(mdqcv3_000 == 0);
    REQUIRE(mdqcv3_001 == 1);
    REQUIRE(mdqcv3_002 == 2);
    REQUIRE(mdqcv3_010 == 3);
    REQUIRE(mdqcv3_011 == 4);
    REQUIRE(mdqcv3_012 == 5);

    int*** ptr = mdqcv3.ptr();
    ptr[0][0][0] = 4;
    ptr[0][0][1] = 5;
    ptr[0][0][2] = 6;
    ptr[0][1][0] = 7;
    ptr[0][1][1] = 8;
    ptr[0][1][2] = 9;
    mdqcv3_000 = ptr[0][0][0];
    mdqcv3_001 = ptr[0][0][1];
    mdqcv3_002 = ptr[0][0][2];
    mdqcv3_010 = ptr[0][1][0];
    mdqcv3_011 = ptr[0][1][1];
    mdqcv3_012 = ptr[0][1][2];
    REQUIRE(mdqcv3[0][0][0] == 4);
    REQUIRE(mdqcv3[0][0][1] == 5);
    REQUIRE(mdqcv3[0][0][2] == 6);
    REQUIRE(mdqcv3[0][1][0] == 7);
    REQUIRE(mdqcv3[0][1][1] == 8);
    REQUIRE(mdqcv3[0][1][2] == 9);
    REQUIRE(mdqcv3_000 == 4);
    REQUIRE(mdqcv3_001 == 5);
    REQUIRE(mdqcv3_002 == 6);
    REQUIRE(mdqcv3_010 == 7);
    REQUIRE(mdqcv3_011 == 8);
    REQUIRE(mdqcv3_012 == 9);

    mdqcv3.deallocate();
    REQUIRE(mdqcv3.size() == 0);
    REQUIRE(mdqcv3.size2() == 0);
    REQUIRE(mdqcv3.size3() == 0);
}

TEST_CASE("darray-6", "[MDQCVector4D]") {
    MDQCVector4D<int> mdqcv4;
    REQUIRE(mdqcv4.size() == 0);
    REQUIRE(mdqcv4.size2() == 0);
    REQUIRE(mdqcv4.size3() == 0);
    REQUIRE(mdqcv4.size4() == 0);

    mdqcv4.allocate(2,3,4,5);
    REQUIRE(mdqcv4.size() == 2);
    REQUIRE(mdqcv4.size2() == 3);
    REQUIRE(mdqcv4.size3() == 4);
    REQUIRE(mdqcv4.size4() == 5);

    mdqcv4.allocate(2,1,1,3);
    REQUIRE(mdqcv4.size() == 2);
    REQUIRE(mdqcv4.size2() == 1);
    REQUIRE(mdqcv4.size3() == 1);
    REQUIRE(mdqcv4.size4() == 3);

    mdqcv4[0][0][0][0] = 0;
    mdqcv4[0][0][0][1] = 1;
    mdqcv4[0][0][0][2] = 2;
    mdqcv4[1][0][0][0] = 3;
    mdqcv4[1][0][0][1] = 4;
    mdqcv4[1][0][0][2] = 5;
    int mdqcv4_0000 = mdqcv4[0][0][0][0];
    int mdqcv4_0001 = mdqcv4[0][0][0][1];
    int mdqcv4_0002 = mdqcv4[0][0][0][2];
    int mdqcv4_1000 = mdqcv4[1][0][0][0];
    int mdqcv4_1001 = mdqcv4[1][0][0][1];
    int mdqcv4_1002 = mdqcv4[1][0][0][2];
    REQUIRE(mdqcv4[0][0][0][0] == 0);
    REQUIRE(mdqcv4[0][0][0][1] == 1);
    REQUIRE(mdqcv4[0][0][0][2] == 2);
    REQUIRE(mdqcv4[1][0][0][0] == 3);
    REQUIRE(mdqcv4[1][0][0][1] == 4);
    REQUIRE(mdqcv4[1][0][0][2] == 5);
    REQUIRE(mdqcv4_0000 == 0);
    REQUIRE(mdqcv4_0001 == 1);
    REQUIRE(mdqcv4_0002 == 2);
    REQUIRE(mdqcv4_1000 == 3);
    REQUIRE(mdqcv4_1001 == 4);
    REQUIRE(mdqcv4_1002 == 5);

    int**** ptr = mdqcv4.ptr();
    ptr[0][0][0][0] = 4;
    ptr[0][0][0][1] = 5;
    ptr[0][0][0][2] = 6;
    ptr[1][0][0][0] = 7;
    ptr[1][0][0][1] = 8;
    ptr[1][0][0][2] = 9;
    mdqcv4_0000 = ptr[0][0][0][0];
    mdqcv4_0001 = ptr[0][0][0][1];
    mdqcv4_0002 = ptr[0][0][0][2];
    mdqcv4_1000 = ptr[1][0][0][0];
    mdqcv4_1001 = ptr[1][0][0][1];
    mdqcv4_1002 = ptr[1][0][0][2];
    REQUIRE(mdqcv4[0][0][0][0] == 4);
    REQUIRE(mdqcv4[0][0][0][1] == 5);
    REQUIRE(mdqcv4[0][0][0][2] == 6);
    REQUIRE(mdqcv4[1][0][0][0] == 7);
    REQUIRE(mdqcv4[1][0][0][1] == 8);
    REQUIRE(mdqcv4[1][0][0][2] == 9);
    REQUIRE(mdqcv4_0000 == 4);
    REQUIRE(mdqcv4_0001 == 5);
    REQUIRE(mdqcv4_0002 == 6);
    REQUIRE(mdqcv4_1000 == 7);
    REQUIRE(mdqcv4_1001 == 8);
    REQUIRE(mdqcv4_1002 == 9);

    mdqcv4.deallocate();
    REQUIRE(mdqcv4.size() == 0);
    REQUIRE(mdqcv4.size2() == 0);
    REQUIRE(mdqcv4.size3() == 0);
    REQUIRE(mdqcv4.size4() == 0);
}
}
