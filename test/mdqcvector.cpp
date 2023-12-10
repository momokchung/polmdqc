#include "catch.hpp"
#include "mdqcvector.h"

using namespace mdqcvector;

TEST_CASE("mdqcvector") {
    Vector1D<int> mdqcv11;
    REQUIRE(mdqcv11.size() == 0);

    mdqcv11.allocate(10);
    REQUIRE(mdqcv11.size() == 10);

    mdqcv11.allocate(5);
    REQUIRE(mdqcv11.size() == 5);

    mdqcv11[0] = 0;
    mdqcv11[1] = 1;
    mdqcv11[2] = 2;
    mdqcv11[3] = 3;
    mdqcv11[4] = 4;
    int mdqcv11_0 = mdqcv11[0];
    int mdqcv11_1 = mdqcv11[1];
    int mdqcv11_2 = mdqcv11[2];
    int mdqcv11_3 = mdqcv11[3];
    int mdqcv11_4 = mdqcv11[4];
    REQUIRE(mdqcv11[0] == 0);
    REQUIRE(mdqcv11[1] == 1);
    REQUIRE(mdqcv11[2] == 2);
    REQUIRE(mdqcv11[3] == 3);
    REQUIRE(mdqcv11[4] == 4);
    REQUIRE(mdqcv11_0 == 0);
    REQUIRE(mdqcv11_1 == 1);
    REQUIRE(mdqcv11_2 == 2);
    REQUIRE(mdqcv11_3 == 3);
    REQUIRE(mdqcv11_4 == 4);

    mdqcv11.deallocate();
    REQUIRE(mdqcv11.size() == 0);

    Vector1D<int> mdqcv12(5);
    mdqcv12[0] = 4;
    mdqcv12[1] = 3;
    mdqcv12[2] = 2;
    mdqcv12[3] = 1;
    mdqcv12[4] = 0;
    REQUIRE(mdqcv12.size() == 5);
    REQUIRE(mdqcv12[0] == 4);
    REQUIRE(mdqcv12[1] == 3);
    REQUIRE(mdqcv12[2] == 2);
    REQUIRE(mdqcv12[3] == 1);
    REQUIRE(mdqcv12[4] == 0);

    Vector2D<double> mdqcv21;
    mdqcv21.allocate(5,6);
    REQUIRE(mdqcv21.size() == 5);
    REQUIRE(mdqcv21.numRows() == 5);
    REQUIRE(mdqcv21.numCols() == 6);

    mdqcv21.allocate(2,3);
    REQUIRE(mdqcv21.size() == 2);
    REQUIRE(mdqcv21.numRows() == 2);
    REQUIRE(mdqcv21.numCols() == 3);
    mdqcv21[0][0] = 0.;
    mdqcv21[0][1] = 1.;
    mdqcv21[0][2] = 2.;
    mdqcv21[1][0] = 3.;
    mdqcv21[1][1] = 4.;
    mdqcv21[1][2] = 5.;

    double mdqcv21_00 = mdqcv21[0][0];
    double mdqcv21_01 = mdqcv21[0][1];
    double mdqcv21_02 = mdqcv21[0][2];
    double mdqcv21_10 = mdqcv21[1][0];
    double mdqcv21_11 = mdqcv21[1][1];
    double mdqcv21_12 = mdqcv21[1][2];

    REQUIRE(mdqcv21_00 == 0.);
    REQUIRE(mdqcv21_01 == 1.);
    REQUIRE(mdqcv21_02 == 2.);
    REQUIRE(mdqcv21_10 == 3.);
    REQUIRE(mdqcv21_11 == 4.);
    REQUIRE(mdqcv21_12 == 5.);

    REQUIRE(mdqcv21[0][0] == 0.);
    REQUIRE(mdqcv21[0][1] == 1.);
    REQUIRE(mdqcv21[0][2] == 2.);
    REQUIRE(mdqcv21[1][0] == 3.);
    REQUIRE(mdqcv21[1][1] == 4.);
    REQUIRE(mdqcv21[1][2] == 5.);

    mdqcv21.deallocate();
    REQUIRE(mdqcv21.size() == 0);
    REQUIRE(mdqcv21.numRows() == 0);
    REQUIRE(mdqcv21.numCols() == 0);

    Vector2D<double> mdqcv22(2,3);
    REQUIRE(mdqcv22.size() == 2);
    REQUIRE(mdqcv22.numRows() == 2);
    REQUIRE(mdqcv22.numCols() == 3);
    mdqcv22[0][0] = 5.;
    mdqcv22[0][1] = 4.;
    mdqcv22[0][2] = 3.;
    mdqcv22[1][0] = 2.;
    mdqcv22[1][1] = 1.;
    mdqcv22[1][2] = 0.;
    REQUIRE(mdqcv22[0][0] == 5.);
    REQUIRE(mdqcv22[0][1] == 4.);
    REQUIRE(mdqcv22[0][2] == 3.);
    REQUIRE(mdqcv22[1][0] == 2.);
    REQUIRE(mdqcv22[1][1] == 1.);
    REQUIRE(mdqcv22[1][2] == 0.);

    Vector3D<int> mdqcv3;
    mdqcv3.allocate(5,6,7);
    REQUIRE(mdqcv3.size() == 5);
    REQUIRE(mdqcv3.sizeX() == 5);
    REQUIRE(mdqcv3.sizeY() == 6);
    REQUIRE(mdqcv3.sizeZ() == 7);

    mdqcv3.allocate(0,0,0);
    REQUIRE(mdqcv3.size() == 0);
    REQUIRE(mdqcv3.sizeX() == 0);
    REQUIRE(mdqcv3.sizeY() == 0);
    REQUIRE(mdqcv3.sizeZ() == 0);

    mdqcv3.allocate(1,2,3);
    REQUIRE(mdqcv3.size() == 1);
    REQUIRE(mdqcv3.sizeX() == 1);
    REQUIRE(mdqcv3.sizeY() == 2);
    REQUIRE(mdqcv3.sizeZ() == 3);
    mdqcv3[0][0][0] = 0;
    mdqcv3[0][0][1] = 1;
    mdqcv3[0][0][2] = 2;
    mdqcv3[0][1][0] = 3;
    mdqcv3[0][1][1] = 4;
    mdqcv3[0][1][2] = 5;

    REQUIRE(mdqcv3[0][0][0] == 0);
    REQUIRE(mdqcv3[0][0][1] == 1);
    REQUIRE(mdqcv3[0][0][2] == 2);
    REQUIRE(mdqcv3[0][1][0] == 3);
    REQUIRE(mdqcv3[0][1][1] == 4);
    REQUIRE(mdqcv3[0][1][2] == 5);

    Vector4D<int> mdqcv4;
    mdqcv4.allocate(2,1,1,3);
    REQUIRE(mdqcv4.size() == 2);
    REQUIRE(mdqcv4.sizeW() == 2);
    REQUIRE(mdqcv4.sizeX() == 1);
    REQUIRE(mdqcv4.sizeY() == 1);
    REQUIRE(mdqcv4.sizeZ() == 3);
    mdqcv4[0][0][0][0] = 0;
    mdqcv4[0][0][0][1] = 1;
    mdqcv4[0][0][0][2] = 2;
    mdqcv4[1][0][0][0] = 3;
    mdqcv4[1][0][0][1] = 4;
    mdqcv4[1][0][0][2] = 5;

    REQUIRE(mdqcv4[0][0][0][0] == 0);
    REQUIRE(mdqcv4[0][0][0][1] == 1);
    REQUIRE(mdqcv4[0][0][0][2] == 2);
    REQUIRE(mdqcv4[1][0][0][0] == 3);
    REQUIRE(mdqcv4[1][0][0][1] == 4);
    REQUIRE(mdqcv4[1][0][0][2] == 5);
}
