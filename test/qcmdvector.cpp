#include "catch.hpp"
#include "qcmdvector.h"

TEST_CASE("qcmdvector") {
    Vector1D<int> qcmdv11;
    REQUIRE(qcmdv11.size() == 0);

    qcmdv11.allocate(10);
    REQUIRE(qcmdv11.size() == 10);

    qcmdv11.allocate(5);
    REQUIRE(qcmdv11.size() == 5);

    qcmdv11[0] = 0;
    qcmdv11[1] = 1;
    qcmdv11[2] = 2;
    qcmdv11[3] = 3;
    qcmdv11[4] = 4;
    int qcmdv11_0 = qcmdv11[0];
    int qcmdv11_1 = qcmdv11[1];
    int qcmdv11_2 = qcmdv11[2];
    int qcmdv11_3 = qcmdv11[3];
    int qcmdv11_4 = qcmdv11[4];
    REQUIRE(qcmdv11[0] == 0);
    REQUIRE(qcmdv11[1] == 1);
    REQUIRE(qcmdv11[2] == 2);
    REQUIRE(qcmdv11[3] == 3);
    REQUIRE(qcmdv11[4] == 4);
    REQUIRE(qcmdv11_0 == 0);
    REQUIRE(qcmdv11_1 == 1);
    REQUIRE(qcmdv11_2 == 2);
    REQUIRE(qcmdv11_3 == 3);
    REQUIRE(qcmdv11_4 == 4);

    qcmdv11.deallocate();
    REQUIRE(qcmdv11.size() == 0);

    Vector1D<int> qcmdv12(5);
    qcmdv12[0] = 4;
    qcmdv12[1] = 3;
    qcmdv12[2] = 2;
    qcmdv12[3] = 1;
    qcmdv12[4] = 0;
    REQUIRE(qcmdv12.size() == 5);
    REQUIRE(qcmdv12[0] == 4);
    REQUIRE(qcmdv12[1] == 3);
    REQUIRE(qcmdv12[2] == 2);
    REQUIRE(qcmdv12[3] == 1);
    REQUIRE(qcmdv12[4] == 0);

    Vector2D<double> qcmdv21;
    qcmdv21.allocate(5,6);
    REQUIRE(qcmdv21.size() == 5);
    REQUIRE(qcmdv21.numRows() == 5);
    REQUIRE(qcmdv21.numCols() == 6);

    qcmdv21.allocate(2,3);
    REQUIRE(qcmdv21.size() == 2);
    REQUIRE(qcmdv21.numRows() == 2);
    REQUIRE(qcmdv21.numCols() == 3);
    qcmdv21[0][0] = 0.;
    qcmdv21[0][1] = 1.;
    qcmdv21[0][2] = 2.;
    qcmdv21[1][0] = 3.;
    qcmdv21[1][1] = 4.;
    qcmdv21[1][2] = 5.;

    double qcmdv21_00 = qcmdv21[0][0];
    double qcmdv21_01 = qcmdv21[0][1];
    double qcmdv21_02 = qcmdv21[0][2];
    double qcmdv21_10 = qcmdv21[1][0];
    double qcmdv21_11 = qcmdv21[1][1];
    double qcmdv21_12 = qcmdv21[1][2];

    REQUIRE(qcmdv21_00 == 0.);
    REQUIRE(qcmdv21_01 == 1.);
    REQUIRE(qcmdv21_02 == 2.);
    REQUIRE(qcmdv21_10 == 3.);
    REQUIRE(qcmdv21_11 == 4.);
    REQUIRE(qcmdv21_12 == 5.);

    REQUIRE(qcmdv21[0][0] == 0.);
    REQUIRE(qcmdv21[0][1] == 1.);
    REQUIRE(qcmdv21[0][2] == 2.);
    REQUIRE(qcmdv21[1][0] == 3.);
    REQUIRE(qcmdv21[1][1] == 4.);
    REQUIRE(qcmdv21[1][2] == 5.);

    qcmdv21.deallocate();
    REQUIRE(qcmdv21.size() == 0);
    REQUIRE(qcmdv21.numRows() == 0);
    REQUIRE(qcmdv21.numCols() == 0);

    Vector2D<double> qcmdv22(2,3);
    REQUIRE(qcmdv22.size() == 2);
    REQUIRE(qcmdv22.numRows() == 2);
    REQUIRE(qcmdv22.numCols() == 3);
    qcmdv22[0][0] = 5.;
    qcmdv22[0][1] = 4.;
    qcmdv22[0][2] = 3.;
    qcmdv22[1][0] = 2.;
    qcmdv22[1][1] = 1.;
    qcmdv22[1][2] = 0.;
    REQUIRE(qcmdv22[0][0] == 5.);
    REQUIRE(qcmdv22[0][1] == 4.);
    REQUIRE(qcmdv22[0][2] == 3.);
    REQUIRE(qcmdv22[1][0] == 2.);
    REQUIRE(qcmdv22[1][1] == 1.);
    REQUIRE(qcmdv22[1][2] == 0.);

    Vector3D<int> qcmdv3;
    qcmdv3.allocate(5,6,7);
    REQUIRE(qcmdv3.size() == 5);
    REQUIRE(qcmdv3.sizeX() == 5);
    REQUIRE(qcmdv3.sizeY() == 6);
    REQUIRE(qcmdv3.sizeZ() == 7);

    qcmdv3.allocate(0,0,0);
    REQUIRE(qcmdv3.size() == 0);
    REQUIRE(qcmdv3.sizeX() == 0);
    REQUIRE(qcmdv3.sizeY() == 0);
    REQUIRE(qcmdv3.sizeZ() == 0);

    qcmdv3.allocate(1,2,3);
    REQUIRE(qcmdv3.size() == 1);
    REQUIRE(qcmdv3.sizeX() == 1);
    REQUIRE(qcmdv3.sizeY() == 2);
    REQUIRE(qcmdv3.sizeZ() == 3);
    qcmdv3[0][0][0] = 0;
    qcmdv3[0][0][1] = 1;
    qcmdv3[0][0][2] = 2;
    qcmdv3[0][1][0] = 3;
    qcmdv3[0][1][1] = 4;
    qcmdv3[0][1][2] = 5;

    REQUIRE(qcmdv3[0][0][0] == 0);
    REQUIRE(qcmdv3[0][0][1] == 1);
    REQUIRE(qcmdv3[0][0][2] == 2);
    REQUIRE(qcmdv3[0][1][0] == 3);
    REQUIRE(qcmdv3[0][1][1] == 4);
    REQUIRE(qcmdv3[0][1][2] == 5);

    Vector4D<int> qcmdv4;
    qcmdv4.allocate(2,1,1,3);
    REQUIRE(qcmdv4.size() == 2);
    REQUIRE(qcmdv4.sizeW() == 2);
    REQUIRE(qcmdv4.sizeX() == 1);
    REQUIRE(qcmdv4.sizeY() == 1);
    REQUIRE(qcmdv4.sizeZ() == 3);
    qcmdv4[0][0][0][0] = 0;
    qcmdv4[0][0][0][1] = 1;
    qcmdv4[0][0][0][2] = 2;
    qcmdv4[1][0][0][0] = 3;
    qcmdv4[1][0][0][1] = 4;
    qcmdv4[1][0][0][2] = 5;

    REQUIRE(qcmdv4[0][0][0][0] == 0);
    REQUIRE(qcmdv4[0][0][0][1] == 1);
    REQUIRE(qcmdv4[0][0][0][2] == 2);
    REQUIRE(qcmdv4[1][0][0][0] == 3);
    REQUIRE(qcmdv4[1][0][0][1] == 4);
    REQUIRE(qcmdv4[1][0][0][2] == 5);
}
